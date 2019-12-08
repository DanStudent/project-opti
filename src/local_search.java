static int[] localSearch(int[] current){
	int numSommet = sommet(current.length);
	int[] best = Arrays.copyOf(current, current.length);
	int bestSolution = getSolution(current);

	int x = Random().nextInt(current.length);

	int y = Random().nextInt(current.length);
	while(x == y){
		y = Random().nextInt(current.length);
	}
	int z = current[x];
	current[x] = current[y];
	current[y] = current[x];

	int newSolution = getSolution(current);
	if(newSolution < bestSolution){
		bestSolution = newSolution;
		return current;
	}
	else{
		return best;
	}
}