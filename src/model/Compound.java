package model;

public class Compound {
	String name;
	
	public Compound(String name){
		this.name = name;
	}

	@Override
	public boolean equals(Object obj) {
		if(obj instanceof Compound){
			return this.name.equals(((Compound) obj).name);
		} else {
			return false;
		}
	}
}
