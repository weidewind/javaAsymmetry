
public class Chromosome {
	private int size;
	private int ori_start;
	private int ori_end;
	private int ter_start;
	private int ter_end;
	private boolean circular;
	private boolean multiori;
	private int configuration;
	
	public Chromosome(int size, int ori_start, int ori_end, int ter_start, int ter_end, boolean circular){
		this.size = size;
		this.ori_start = ori_start;
		this.ori_end = ori_end;
		this.ter_start = ter_start;
		this.ter_end = ter_end;
		this.circular = circular;
		this.setConfiguration();
	}
	
	public Chromosome(int size, int ori_start, int ori_end, int ter, boolean circular, boolean multiori){
		this.size = size;
		this.ori_start = ori_start;
		this.ori_end = ori_end;
		this.ter_start = ter;
		this.ter_end = ter;
		this.circular = circular;
		this.multiori = multiori;
		this.setConfiguration();
	}
	


	public void setOri_start(int ori_start) {
		this.ori_start = ori_start;
	}

	public void setOri_end(int ori_end) {
		this.ori_end = ori_end;
	}

	public void setTer_start(int ter_start) {
		this.ter_start = ter_start;
	}

	public void setTer_end(int ter_end) {
		this.ter_end = ter_end;
	}
	
	public void setConfiguration(int configuration) {
		this.configuration = configuration;
	}


	public int getSize() {
		return size;
	}

	public int getOri_start() {
		return ori_start;
	}

	public int getOri_end() {
		return ori_end;
	}

	public int getTer_start() {
		return ter_start;
	}
	
	public int getTer_end() {
		return ter_end;
	}
	
	
	public int getConfiguration(){ 
		return configuration;
	}

	public boolean isCircular() {
		return circular;
	}
	
	public boolean isMultiOri() {
		return multiori;
	}
	
/** Expands ori region by _length_ nucleotides in both directions
 * 	
 * @param length
 */
	public void expandOri(int length){
		if (this.configuration == 1){
			if (this.ori_start >= length){
				this.setOri_start(this.ori_start - length);
				this.setOri_end(this.ori_end + length);
			}
			else {
				this.setOri_start(this.size + this.ori_start - length);
				this.setOri_end(this.ori_end + length);
			}
		}
		if (this.configuration == 2 | this.configuration == 4){
			this.setOri_start(this.ori_start - length);
			this.setOri_end(this.ori_end + length);
		}
		if (this.configuration == 3){
			if (this.ori_end + length <= this.size){
				this.setOri_start(this.ori_start - length);
				this.setOri_end(this.ori_end + length);
			}
			else {
				this.setOri_start(this.ori_start - length);
				this.setOri_end(length + this.ori_end - this.size);
			}
		}
		this.setConfiguration();
	}
	
	public void expandTer(int length){
		if (this.configuration == 3){
			if (this.ter_start >= length){
				this.setTer_start(this.ter_start - length);
				this.setTer_end(this.ter_end + length);
			}
			else {
				this.setTer_start(this.size + this.ter_start - length);
				this.setTer_end(this.ter_end + length);
			}
		}
		if (this.configuration == 2 | this.configuration == 4){
			this.setTer_start(this.ter_start - length);
			this.setTer_end(this.ter_end + length);
		}
		if (this.configuration == 1){
			if (this.ter_end + length <= this.size){
				this.setTer_start(this.ter_start - length);
				this.setTer_end(this.ter_end + length);
			}
			else {
				this.setTer_start(this.ter_start - length);
				this.setTer_end(length + this.ter_end - this.size);
			}
		}
		this.setConfiguration();
	}
	
	public void setConfiguration(){
		
		if (this.ori_end < this.ter_start){
			if (this.ori_start < this.ter_end){
				this.configuration = 1;
			}
			else {
				if(this.ter_start <= this.ter_end){
					this.configuration = 2;
				}
				else this.configuration  = 4;
			}
		}
		else this.configuration = 3;
	}
	
}