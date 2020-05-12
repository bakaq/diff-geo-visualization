#include <SDL2/SDL.h>

#include <array>
#include <math.h>
#include <stdio.h>

const int SCREEN_W = 640;
const int SCREEN_H = 480;

const int SCALE = 100;

const double PI = 3.141592;

// Records the coordinates of the vector and it's components
// (basically the components in the tangent bundle I guess)
struct Vec{
	std::array<double, 2> pos = {0,0};
	std::array<double, 2> vel = {0,0};
};


// Define metric
typedef std::array<std::array<double, 2>, 2> Matrix2;

Matrix2 metric(std::array<double,2> coord){
	Matrix2 m = {std::array<double,2>{1,0},
				 std::array<double,2>{0,pow(coord[0], 2)}};
	return m;
}

Matrix2 inv_metric(std::array<double,2> coord){
	// TODO: Non diagonal metrics
	Matrix2 m = metric(coord);
	Matrix2 inv_m = {std::array<double,2>{1/m[0][0], 0},
					 std::array<double,2>{0,1/m[1][1]}};
	return inv_m;
}

// Get Christoffel Symbols
typedef std::array<Matrix2, 2> Matrix3;

Matrix3 gamma(std::array<double, 2> coord){
	double h = 0.000001;
	Matrix3 Gamma;

	for(int i = 0; i < 2; ++i){
	 for(int j = 0; j < 2; ++j){
	  for(int k = 0; k < 2; ++k){
		//printf("Calculating gamma: i,j,k = %d, %d, %d\n", i, j, k);
		if(k > j){
			// Symetrization
			Gamma[i][j][k] = Gamma[i][k][j];		
		}else{
			Gamma[i][j][k] = 0;
			
			for(int l = 0; l < 2; ++l){
				//printf("l = %d\n", l);
				std::array<double, 2> delta1 = coord, delta2 = coord, delta3 = coord;
				
				delta1[k] += h;
				double term1 = (metric(delta1)[l][j] - metric(coord)[l][j])/h;

				delta2[j] += h;
				double term2 = (metric(delta2)[l][k] - metric(coord)[l][k])/h;

				delta3[l] += h;
				double term3 = (metric(delta3)[j][k] - metric(coord)[j][k])/h;
				
				Gamma[i][j][k] += inv_metric(coord)[i][l]*(term1 + term2 - term3);
			}
			Gamma[i][j][k] = (0.5)*Gamma[i][j][k];
		}
	  }
	 }
	}	

	return Gamma;
}


int main(){
	
	// Windows n shit
	SDL_Window *window = NULL;
	SDL_Renderer *renderer = NULL;
	//SDL_Texture *texture = NULL;

	// Init SDL
	SDL_Init(SDL_INIT_VIDEO);

	// Initialize windows n shit
	window = SDL_CreateWindow("Differential Geometry Visualization",
			SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
			SCREEN_W, SCREEN_H, SDL_WINDOW_SHOWN);

	renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
	SDL_SetRenderDrawColor(renderer, 0xFF,0xFF,0xFF,0xFF);

	// Initialize vector
	Vec v;
	v.pos = {1,0};
	v.vel = {0,0.5};

	// Main loop
	bool running = true;
	SDL_Event e;
	while(running){
	
		// Events
		while(SDL_PollEvent(&e) != 0){
			if(e.type == SDL_QUIT){
				running = false;
			}
		}
		
		// Translate vector
		std::array<double, 2> dx = {0,0.001};
		std::array<double, 2> dv;
		Matrix3 G = gamma(v.pos);
		for(int i = 0; i < 2; ++i){
			dv[i] = 0;
			for(int j = 0; j < 2; ++j){
				for(int k = 0; k < 2; ++k){
					dv[i] -= G[i][j][k]*v.vel[j]*dx[k];
				}
			}
		}
		
		for(int i = 0; i < 2; ++i){
			v.pos[i] += dx[i];
			v.vel[i] += dv[i];
		}

		// Clear renderer 
		SDL_SetRenderDrawColor(renderer, 0xFF,0xFF,0xFF,0xFF);
		SDL_RenderClear(renderer);
		
		// Draw grid
		SDL_SetRenderDrawColor(renderer, 0x00,0x00,0x00,0xFF);
		SDL_RenderDrawLine(renderer,
				0, SCREEN_H/2,SCREEN_W, SCREEN_H/2);
		SDL_RenderDrawLine(renderer,
				SCREEN_W/2, 0, SCREEN_W/2, SCREEN_H);
		
		// Draw unit circle (there must be a better way to do this)
		SDL_SetRenderDrawColor(renderer, 0x00,0x00,0xFF,0xFF);
		for(double theta = 0; theta < 2*PI; theta += 1/((double)SCALE)){
			SDL_RenderDrawPoint(renderer,
					SCREEN_W/2 + (int)(cos(theta)*SCALE),
					SCREEN_H/2 + (int)(sin(theta)*SCALE));
		}
		
		// Draw vector
		SDL_SetRenderDrawColor(renderer, 0xFF,0x00,0x00,0xFF);
		SDL_RenderDrawLine(renderer,
			SCREEN_W/2 + (int)(v.pos[0]*cos(v.pos[1])*SCALE),
			SCREEN_H/2 - (int)(v.pos[0]*sin(v.pos[1])*SCALE),
			SCREEN_W/2 + (int)((v.pos[0]*cos(v.pos[1])-v.pos[0]*v.vel[1]*sin(v.pos[1])+v.vel[0]*cos(v.pos[1]))*SCALE),
			SCREEN_H/2 - (int)((v.pos[0]*sin(v.pos[1])+v.pos[0]*v.vel[1]*cos(v.pos[1])+v.vel[0]*sin(v.pos[1]))*SCALE));

		SDL_RenderPresent(renderer);
		
		// Frame cap
		SDL_Delay((int)(1000/600));
	}

	// Free stuff
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	SDL_Quit();
}
