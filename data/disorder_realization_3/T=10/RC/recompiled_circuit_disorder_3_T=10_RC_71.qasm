OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.91575032) q[0];
sx q[0];
rz(-3.1103818) q[0];
sx q[0];
rz(-2.6565235) q[0];
rz(0.78753161) q[1];
sx q[1];
rz(-1.0163611) q[1];
sx q[1];
rz(5.8689868) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.480455) q[0];
sx q[0];
rz(-1.2167131) q[0];
sx q[0];
rz(-2.7829091) q[0];
x q[1];
rz(-1.4235052) q[2];
sx q[2];
rz(-2.5669332) q[2];
sx q[2];
rz(1.6342083) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.011885) q[1];
sx q[1];
rz(-1.1457448) q[1];
sx q[1];
rz(-0.070406291) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6730509) q[3];
sx q[3];
rz(-1.8342606) q[3];
sx q[3];
rz(-1.9139569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9238613) q[2];
sx q[2];
rz(-1.2552746) q[2];
sx q[2];
rz(3.1100173) q[2];
rz(1.2565553) q[3];
sx q[3];
rz(-2.6387408) q[3];
sx q[3];
rz(0.52662915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8388222) q[0];
sx q[0];
rz(-1.6844203) q[0];
sx q[0];
rz(0.1698499) q[0];
rz(-0.70392144) q[1];
sx q[1];
rz(-2.070065) q[1];
sx q[1];
rz(-2.6020715) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9233421) q[0];
sx q[0];
rz(-1.6174416) q[0];
sx q[0];
rz(1.1211066) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1621446) q[2];
sx q[2];
rz(-2.2505629) q[2];
sx q[2];
rz(-2.0603927) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9425548) q[1];
sx q[1];
rz(-2.6032762) q[1];
sx q[1];
rz(1.1953137) q[1];
rz(-pi) q[2];
rz(2.1971365) q[3];
sx q[3];
rz(-0.8144905) q[3];
sx q[3];
rz(-2.8163547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4743621) q[2];
sx q[2];
rz(-1.2381866) q[2];
sx q[2];
rz(-0.24307069) q[2];
rz(2.4754751) q[3];
sx q[3];
rz(-0.56454286) q[3];
sx q[3];
rz(-1.8977785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3617525) q[0];
sx q[0];
rz(-3.0267974) q[0];
sx q[0];
rz(0.4483805) q[0];
rz(1.7547296) q[1];
sx q[1];
rz(-1.153839) q[1];
sx q[1];
rz(0.2562491) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8921593) q[0];
sx q[0];
rz(-1.8349577) q[0];
sx q[0];
rz(0.63412068) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44552866) q[2];
sx q[2];
rz(-2.0085213) q[2];
sx q[2];
rz(0.23756269) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8505893) q[1];
sx q[1];
rz(-0.91445078) q[1];
sx q[1];
rz(-2.7235051) q[1];
rz(-2.9647397) q[3];
sx q[3];
rz(-1.5336509) q[3];
sx q[3];
rz(2.0166486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3391352) q[2];
sx q[2];
rz(-0.70636237) q[2];
sx q[2];
rz(0.57470542) q[2];
rz(-1.7859219) q[3];
sx q[3];
rz(-1.1698497) q[3];
sx q[3];
rz(-1.908196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.213585) q[0];
sx q[0];
rz(-1.7127697) q[0];
sx q[0];
rz(0.25948778) q[0];
rz(1.9909987) q[1];
sx q[1];
rz(-1.3535627) q[1];
sx q[1];
rz(2.4096699) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4703092) q[0];
sx q[0];
rz(-1.0571612) q[0];
sx q[0];
rz(1.0259823) q[0];
rz(-pi) q[1];
rz(3.0078997) q[2];
sx q[2];
rz(-2.4198654) q[2];
sx q[2];
rz(2.0267817) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.96958292) q[1];
sx q[1];
rz(-1.3130377) q[1];
sx q[1];
rz(-1.7715363) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51283522) q[3];
sx q[3];
rz(-0.25179112) q[3];
sx q[3];
rz(2.4076715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4449473) q[2];
sx q[2];
rz(-1.8680633) q[2];
sx q[2];
rz(2.990492) q[2];
rz(-2.5949196) q[3];
sx q[3];
rz(-1.0497382) q[3];
sx q[3];
rz(-2.7643519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5916409) q[0];
sx q[0];
rz(-2.0786091) q[0];
sx q[0];
rz(-1.8792101) q[0];
rz(1.6732015) q[1];
sx q[1];
rz(-0.60931283) q[1];
sx q[1];
rz(2.343822) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7933465) q[0];
sx q[0];
rz(-3.0214546) q[0];
sx q[0];
rz(-1.5500463) q[0];
rz(-pi) q[1];
rz(1.7902137) q[2];
sx q[2];
rz(-1.2741538) q[2];
sx q[2];
rz(0.23362939) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.080575374) q[1];
sx q[1];
rz(-1.4772381) q[1];
sx q[1];
rz(0.94228014) q[1];
rz(-pi) q[2];
rz(-2.6990715) q[3];
sx q[3];
rz(-1.8707152) q[3];
sx q[3];
rz(-2.4181441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.795934) q[2];
sx q[2];
rz(-0.63085932) q[2];
sx q[2];
rz(-2.8395555) q[2];
rz(-1.1473514) q[3];
sx q[3];
rz(-1.6796422) q[3];
sx q[3];
rz(0.54774493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40925947) q[0];
sx q[0];
rz(-3.0497666) q[0];
sx q[0];
rz(-1.1557895) q[0];
rz(-1.0844768) q[1];
sx q[1];
rz(-0.98025727) q[1];
sx q[1];
rz(0.070080431) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87175831) q[0];
sx q[0];
rz(-1.7239128) q[0];
sx q[0];
rz(1.397875) q[0];
x q[1];
rz(-2.3577865) q[2];
sx q[2];
rz(-1.98181) q[2];
sx q[2];
rz(-2.59936) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8225704) q[1];
sx q[1];
rz(-0.59148568) q[1];
sx q[1];
rz(0.2925847) q[1];
rz(-pi) q[2];
rz(3.1061884) q[3];
sx q[3];
rz(-1.6453214) q[3];
sx q[3];
rz(0.41985928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8391116) q[2];
sx q[2];
rz(-1.9589067) q[2];
sx q[2];
rz(-0.62136674) q[2];
rz(1.7012043) q[3];
sx q[3];
rz(-2.6337603) q[3];
sx q[3];
rz(2.5682209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0550585) q[0];
sx q[0];
rz(-1.4165514) q[0];
sx q[0];
rz(-2.281718) q[0];
rz(-1.2043918) q[1];
sx q[1];
rz(-2.269373) q[1];
sx q[1];
rz(-3.133657) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36695776) q[0];
sx q[0];
rz(-1.8263706) q[0];
sx q[0];
rz(0.5704244) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0646348) q[2];
sx q[2];
rz(-2.3828265) q[2];
sx q[2];
rz(0.90492349) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4850033) q[1];
sx q[1];
rz(-2.2409229) q[1];
sx q[1];
rz(-2.2336002) q[1];
rz(-pi) q[2];
rz(-0.31422024) q[3];
sx q[3];
rz(-1.1952956) q[3];
sx q[3];
rz(0.34208959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4456711) q[2];
sx q[2];
rz(-1.9182703) q[2];
sx q[2];
rz(2.725214) q[2];
rz(-1.3683866) q[3];
sx q[3];
rz(-1.8442644) q[3];
sx q[3];
rz(-0.90464512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1798582) q[0];
sx q[0];
rz(-0.27357736) q[0];
sx q[0];
rz(0.36488786) q[0];
rz(0.94003135) q[1];
sx q[1];
rz(-0.53661984) q[1];
sx q[1];
rz(-1.6392802) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68740326) q[0];
sx q[0];
rz(-0.25035509) q[0];
sx q[0];
rz(-3.1012015) q[0];
rz(-pi) q[1];
rz(2.084923) q[2];
sx q[2];
rz(-0.80386111) q[2];
sx q[2];
rz(0.67509292) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.18174905) q[1];
sx q[1];
rz(-1.0105003) q[1];
sx q[1];
rz(0.77397857) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3650465) q[3];
sx q[3];
rz(-1.4806517) q[3];
sx q[3];
rz(-0.84116018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6442948) q[2];
sx q[2];
rz(-2.6361894) q[2];
sx q[2];
rz(-2.0765182) q[2];
rz(2.8403357) q[3];
sx q[3];
rz(-1.732429) q[3];
sx q[3];
rz(1.6206954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.168468) q[0];
sx q[0];
rz(-0.081806101) q[0];
sx q[0];
rz(2.705943) q[0];
rz(1.3849974) q[1];
sx q[1];
rz(-2.6711617) q[1];
sx q[1];
rz(-0.41697821) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0957085) q[0];
sx q[0];
rz(-2.4864712) q[0];
sx q[0];
rz(1.6326293) q[0];
rz(-0.66395335) q[2];
sx q[2];
rz(-2.0400527) q[2];
sx q[2];
rz(-2.9479153) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4638085) q[1];
sx q[1];
rz(-0.52075547) q[1];
sx q[1];
rz(-1.2893454) q[1];
rz(-pi) q[2];
rz(3.0005089) q[3];
sx q[3];
rz(-1.699563) q[3];
sx q[3];
rz(-1.5878549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.10432648) q[2];
sx q[2];
rz(-1.4929079) q[2];
sx q[2];
rz(2.9157675) q[2];
rz(-0.2078235) q[3];
sx q[3];
rz(-0.72312975) q[3];
sx q[3];
rz(-0.58661714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.726783) q[0];
sx q[0];
rz(-0.2668969) q[0];
sx q[0];
rz(1.6171932) q[0];
rz(-2.1879451) q[1];
sx q[1];
rz(-1.2748268) q[1];
sx q[1];
rz(-1.3226002) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029862558) q[0];
sx q[0];
rz(-2.7435281) q[0];
sx q[0];
rz(-1.8882621) q[0];
rz(-pi) q[1];
x q[1];
rz(0.25090353) q[2];
sx q[2];
rz(-1.13399) q[2];
sx q[2];
rz(1.3383588) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9433371) q[1];
sx q[1];
rz(-2.2554734) q[1];
sx q[1];
rz(2.4556922) q[1];
x q[2];
rz(-0.15354746) q[3];
sx q[3];
rz(-2.2205177) q[3];
sx q[3];
rz(0.49680199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3502729) q[2];
sx q[2];
rz(-2.0952756) q[2];
sx q[2];
rz(-2.1255169) q[2];
rz(-1.2223876) q[3];
sx q[3];
rz(-0.17761579) q[3];
sx q[3];
rz(0.55541742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14324698) q[0];
sx q[0];
rz(-2.2194942) q[0];
sx q[0];
rz(1.0794328) q[0];
rz(1.7779508) q[1];
sx q[1];
rz(-1.2294055) q[1];
sx q[1];
rz(1.3235863) q[1];
rz(-0.50921847) q[2];
sx q[2];
rz(-1.5621395) q[2];
sx q[2];
rz(-0.12315673) q[2];
rz(1.3397459) q[3];
sx q[3];
rz(-1.4828724) q[3];
sx q[3];
rz(2.8433269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
