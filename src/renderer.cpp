#include <limits>
#include <algorithm>
#include "renderer.hpp"

Renderer::Renderer( const fractal_land& land, const pheronome& phen, 
                    const position_t& pos_nest, const position_t& pos_food,
                    const AntPopulation& ants )
    :   m_ref_land( land ),
        m_land( nullptr ),
        m_ref_phen( phen ),
        m_pos_nest( pos_nest ),
        m_pos_food( pos_food ),
        m_ref_ants( ants )
{
}

Renderer::~Renderer() {
    if ( m_land != nullptr )
        SDL_DestroyTexture( m_land );
}

void Renderer::display( Window& win, std::size_t const& compteur )
{
    SDL_Renderer* renderer = SDL_GetRenderer( win.get() );
    
    // Création texture paysage (inchangé)
    if ( m_land == nullptr ) {
        SDL_Surface* temp_surface = SDL_CreateRGBSurface(0, m_ref_land.dimensions(), m_ref_land.dimensions(), 32,
                                                          0x00FF0000, 0x0000FF00, 0x000000FF, 0xFF000000);
        double min_height{std::numeric_limits<double>::max()}, max_height{std::numeric_limits<double>::lowest()};
        for ( fractal_land::dim_t i = 0; i < m_ref_land.dimensions( ); ++i )
            for ( fractal_land::dim_t j = 0; j < m_ref_land.dimensions( ); ++j ) {
                min_height = std::min( min_height, m_ref_land( i, j ) );
                max_height = std::max( max_height, m_ref_land( i, j ) );
            }
        for ( fractal_land::dim_t i = 0; i < m_ref_land.dimensions( ); ++i )
            for ( fractal_land::dim_t j = 0; j < m_ref_land.dimensions( ); ++j ) {
                double c = 255. * ( m_ref_land( i, j ) - min_height ) / ( max_height - min_height );
                Uint32* pixel = (Uint32*) ((Uint8*)temp_surface->pixels + j * temp_surface->pitch + i * sizeof(Uint32));
                *pixel = SDL_MapRGBA( temp_surface->format, static_cast<Uint8>(c), static_cast<Uint8>(c), static_cast<Uint8>(c), 255 );
            }
        m_land = SDL_CreateTextureFromSurface( renderer, temp_surface );
        SDL_FreeSurface( temp_surface );
    }
    
    SDL_SetRenderDrawColor( renderer, 0, 0, 0, 255 );
    SDL_RenderClear( renderer );
    
    SDL_Rect dest_rect1{0, 0, static_cast<int>(m_ref_land.dimensions()), static_cast<int>(m_ref_land.dimensions())};
    SDL_RenderCopy( renderer, m_land, nullptr, &dest_rect1 );
    SDL_Rect dest_rect2{static_cast<int>(m_ref_land.dimensions()) + 10, 0, static_cast<int>(m_ref_land.dimensions()), static_cast<int>(m_ref_land.dimensions())};
    SDL_RenderCopy( renderer, m_land, nullptr, &dest_rect2 );
    
    SDL_SetRenderDrawBlendMode( renderer, SDL_BLENDMODE_BLEND );
    
    // --- MODIFICATION POUR ANT_POPULATION ---
    // On boucle sur la taille des vecteurs
    size_t nb_ants = m_ref_ants.m_x.size();
    win.set_pen( 0, 255, 255 );
    for ( size_t i = 0; i < nb_ants; ++i ) {
        // Accès direct aux tableaux (SoA)
        win.pset( m_ref_ants.m_x[i], m_ref_ants.m_y[i] );
    }
    
    // Affichage phéromones (inchangé)
    for ( fractal_land::dim_t i = 0; i < m_ref_land.dimensions( ); ++i )
        for ( fractal_land::dim_t j = 0; j < m_ref_land.dimensions( ); ++j ) {
            double r = std::min( 1., (double)m_ref_phen( i, j )[0] );
            double g = std::min( 1., (double)m_ref_phen( i, j )[1] );
            if ( r > 0.01 || g > 0.01 ) {
                win.set_pen( static_cast<Uint8>( r * 255 ), static_cast<Uint8>( g * 255 ), 0 );
                win.pset( static_cast<int>( i + m_ref_land.dimensions( ) + 10 ), static_cast<int>( j ) );
            }
        }
    
    // Courbe (inchangé)
    m_curve.push_back(compteur);
    if ( m_curve.size( ) > 1 ) {
        int sz_win = win.size( ).first;
        int ydec = win.size( ).second - 1;
        double max_curve_val = *std::max_element( m_curve.begin(), m_curve.end() );
        double h_max_val = 256. / std::max( max_curve_val, 1.);
        double step      = double(sz_win) / (double)( m_curve.size( ) );
        SDL_SetRenderDrawColor( renderer, 255, 255, 127, 255 );
        for ( std::size_t i = 0; i < m_curve.size( ) - 1; i++ ) {
            int x1 = static_cast<int>( i * step );
            int y1 = static_cast<int>( ydec - m_curve[i] * h_max_val );
            int x2 = static_cast<int>( ( i + 1 ) * step );
            int y2 = static_cast<int>( ydec - m_curve[i + 1] * h_max_val );
            SDL_RenderDrawLine( renderer, x1, y1, x2, y2 );
        }
    }
    SDL_RenderPresent( renderer );
}