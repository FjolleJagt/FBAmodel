package model;

public class Compound {
	public String name;
	public boolean isInMassBalance;

	@Override
	public boolean equals(Object obj) {
		if(obj instanceof Compound){
			return this.name.equals(((Compound) obj).name);
		} else {
			return false;
		}
	}
}
