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
rz(2.7273942) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038923351) q[0];
sx q[0];
rz(-1.9063213) q[0];
sx q[0];
rz(1.9468007) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0010927) q[2];
sx q[2];
rz(-1.6506519) q[2];
sx q[2];
rz(-2.9542838) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.011885) q[1];
sx q[1];
rz(-1.1457448) q[1];
sx q[1];
rz(3.0711864) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6730509) q[3];
sx q[3];
rz(-1.3073321) q[3];
sx q[3];
rz(1.2276358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9238613) q[2];
sx q[2];
rz(-1.2552746) q[2];
sx q[2];
rz(-0.031575354) q[2];
rz(-1.8850373) q[3];
sx q[3];
rz(-0.50285181) q[3];
sx q[3];
rz(-0.52662915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3027705) q[0];
sx q[0];
rz(-1.6844203) q[0];
sx q[0];
rz(-0.1698499) q[0];
rz(2.4376712) q[1];
sx q[1];
rz(-1.0715276) q[1];
sx q[1];
rz(-0.53952113) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9233421) q[0];
sx q[0];
rz(-1.5241511) q[0];
sx q[0];
rz(-1.1211066) q[0];
rz(-pi) q[1];
rz(0.72210724) q[2];
sx q[2];
rz(-1.8849843) q[2];
sx q[2];
rz(-2.9177641) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5121465) q[1];
sx q[1];
rz(-2.0680032) q[1];
sx q[1];
rz(2.9260103) q[1];
rz(-pi) q[2];
rz(2.1971365) q[3];
sx q[3];
rz(-0.8144905) q[3];
sx q[3];
rz(-2.8163547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4743621) q[2];
sx q[2];
rz(-1.903406) q[2];
sx q[2];
rz(2.898522) q[2];
rz(-0.66611755) q[3];
sx q[3];
rz(-2.5770498) q[3];
sx q[3];
rz(-1.2438141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77984017) q[0];
sx q[0];
rz(-3.0267974) q[0];
sx q[0];
rz(0.4483805) q[0];
rz(1.7547296) q[1];
sx q[1];
rz(-1.9877537) q[1];
sx q[1];
rz(2.8853436) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6305144) q[0];
sx q[0];
rz(-0.96195463) q[0];
sx q[0];
rz(1.2468673) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82661144) q[2];
sx q[2];
rz(-0.6140784) q[2];
sx q[2];
rz(-0.6073063) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8505893) q[1];
sx q[1];
rz(-0.91445078) q[1];
sx q[1];
rz(2.7235051) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9647397) q[3];
sx q[3];
rz(-1.6079418) q[3];
sx q[3];
rz(1.1249441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3391352) q[2];
sx q[2];
rz(-0.70636237) q[2];
sx q[2];
rz(-2.5668872) q[2];
rz(1.7859219) q[3];
sx q[3];
rz(-1.1698497) q[3];
sx q[3];
rz(-1.2333966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.213585) q[0];
sx q[0];
rz(-1.7127697) q[0];
sx q[0];
rz(2.8821049) q[0];
rz(1.9909987) q[1];
sx q[1];
rz(-1.78803) q[1];
sx q[1];
rz(-2.4096699) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67128348) q[0];
sx q[0];
rz(-2.0844315) q[0];
sx q[0];
rz(-1.0259823) q[0];
rz(-pi) q[1];
rz(2.4243083) q[2];
sx q[2];
rz(-1.4826164) q[2];
sx q[2];
rz(-0.35536534) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6432453) q[1];
sx q[1];
rz(-2.8162662) q[1];
sx q[1];
rz(-0.64756067) q[1];
rz(-pi) q[2];
rz(-1.4452403) q[3];
sx q[3];
rz(-1.7896277) q[3];
sx q[3];
rz(0.20727508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6966454) q[2];
sx q[2];
rz(-1.8680633) q[2];
sx q[2];
rz(0.15110061) q[2];
rz(2.5949196) q[3];
sx q[3];
rz(-2.0918545) q[3];
sx q[3];
rz(0.37724075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54995173) q[0];
sx q[0];
rz(-1.0629835) q[0];
sx q[0];
rz(-1.2623825) q[0];
rz(-1.6732015) q[1];
sx q[1];
rz(-2.5322798) q[1];
sx q[1];
rz(2.343822) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3691468) q[0];
sx q[0];
rz(-1.6909084) q[0];
sx q[0];
rz(-3.1390879) q[0];
rz(-2.8380978) q[2];
sx q[2];
rz(-1.7804838) q[2];
sx q[2];
rz(1.7393302) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5835727) q[1];
sx q[1];
rz(-0.94545525) q[1];
sx q[1];
rz(0.11548345) q[1];
rz(-1.9005152) q[3];
sx q[3];
rz(-1.1493249) q[3];
sx q[3];
rz(0.7082522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.795934) q[2];
sx q[2];
rz(-2.5107333) q[2];
sx q[2];
rz(0.30203715) q[2];
rz(-1.9942412) q[3];
sx q[3];
rz(-1.6796422) q[3];
sx q[3];
rz(-0.54774493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7323332) q[0];
sx q[0];
rz(-3.0497666) q[0];
sx q[0];
rz(1.9858032) q[0];
rz(-2.0571158) q[1];
sx q[1];
rz(-0.98025727) q[1];
sx q[1];
rz(-0.070080431) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018774059) q[0];
sx q[0];
rz(-0.2304603) q[0];
sx q[0];
rz(-2.3019058) q[0];
x q[1];
rz(0.55307936) q[2];
sx q[2];
rz(-0.86420176) q[2];
sx q[2];
rz(1.4097708) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.97092123) q[1];
sx q[1];
rz(-1.0075924) q[1];
sx q[1];
rz(1.3794273) q[1];
x q[2];
rz(-3.1061884) q[3];
sx q[3];
rz(-1.6453214) q[3];
sx q[3];
rz(-0.41985928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.30248102) q[2];
sx q[2];
rz(-1.9589067) q[2];
sx q[2];
rz(-2.5202259) q[2];
rz(-1.7012043) q[3];
sx q[3];
rz(-2.6337603) q[3];
sx q[3];
rz(-2.5682209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086534111) q[0];
sx q[0];
rz(-1.4165514) q[0];
sx q[0];
rz(0.85987464) q[0];
rz(1.9372008) q[1];
sx q[1];
rz(-2.269373) q[1];
sx q[1];
rz(0.0079356114) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3646274) q[0];
sx q[0];
rz(-1.0210751) q[0];
sx q[0];
rz(-1.2697898) q[0];
x q[1];
rz(2.7107312) q[2];
sx q[2];
rz(-2.2164946) q[2];
sx q[2];
rz(-0.25260392) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4850033) q[1];
sx q[1];
rz(-2.2409229) q[1];
sx q[1];
rz(2.2336002) q[1];
rz(-0.90585917) q[3];
sx q[3];
rz(-2.6568036) q[3];
sx q[3];
rz(-1.0672027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4456711) q[2];
sx q[2];
rz(-1.9182703) q[2];
sx q[2];
rz(2.725214) q[2];
rz(-1.773206) q[3];
sx q[3];
rz(-1.2973283) q[3];
sx q[3];
rz(2.2369475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96173441) q[0];
sx q[0];
rz(-0.27357736) q[0];
sx q[0];
rz(-0.36488786) q[0];
rz(2.2015613) q[1];
sx q[1];
rz(-0.53661984) q[1];
sx q[1];
rz(-1.5023124) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2973328) q[0];
sx q[0];
rz(-1.5808006) q[0];
sx q[0];
rz(-0.25015932) q[0];
rz(-pi) q[1];
rz(0.83606007) q[2];
sx q[2];
rz(-1.2088472) q[2];
sx q[2];
rz(-1.8723633) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8870526) q[1];
sx q[1];
rz(-2.2215543) q[1];
sx q[1];
rz(-0.73144967) q[1];
x q[2];
rz(-3.0495166) q[3];
sx q[3];
rz(-1.7756988) q[3];
sx q[3];
rz(-0.74842194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6442948) q[2];
sx q[2];
rz(-2.6361894) q[2];
sx q[2];
rz(-2.0765182) q[2];
rz(-2.8403357) q[3];
sx q[3];
rz(-1.4091636) q[3];
sx q[3];
rz(-1.5208972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.168468) q[0];
sx q[0];
rz(-3.0597866) q[0];
sx q[0];
rz(-0.43564963) q[0];
rz(1.3849974) q[1];
sx q[1];
rz(-0.47043097) q[1];
sx q[1];
rz(-2.7246144) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1736261) q[0];
sx q[0];
rz(-2.2244503) q[0];
sx q[0];
rz(0.047441479) q[0];
rz(2.4531104) q[2];
sx q[2];
rz(-2.3496029) q[2];
sx q[2];
rz(2.288523) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3560564) q[1];
sx q[1];
rz(-1.0724663) q[1];
sx q[1];
rz(0.15798012) q[1];
x q[2];
rz(1.4407518) q[3];
sx q[3];
rz(-1.4308883) q[3];
sx q[3];
rz(-0.035294447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.10432648) q[2];
sx q[2];
rz(-1.4929079) q[2];
sx q[2];
rz(0.22582516) q[2];
rz(-0.2078235) q[3];
sx q[3];
rz(-2.4184629) q[3];
sx q[3];
rz(0.58661714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41480961) q[0];
sx q[0];
rz(-2.8746958) q[0];
sx q[0];
rz(-1.5243994) q[0];
rz(2.1879451) q[1];
sx q[1];
rz(-1.2748268) q[1];
sx q[1];
rz(1.3226002) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7693217) q[0];
sx q[0];
rz(-1.1936545) q[0];
sx q[0];
rz(-3.0110714) q[0];
x q[1];
rz(-1.0820504) q[2];
sx q[2];
rz(-2.6419123) q[2];
sx q[2];
rz(-2.3479455) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0364914) q[1];
sx q[1];
rz(-2.0836012) q[1];
sx q[1];
rz(0.75863104) q[1];
x q[2];
rz(0.15354746) q[3];
sx q[3];
rz(-2.2205177) q[3];
sx q[3];
rz(2.6447907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7913197) q[2];
sx q[2];
rz(-1.046317) q[2];
sx q[2];
rz(1.0160758) q[2];
rz(-1.919205) q[3];
sx q[3];
rz(-2.9639769) q[3];
sx q[3];
rz(0.55541742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14324698) q[0];
sx q[0];
rz(-0.9220985) q[0];
sx q[0];
rz(-2.0621598) q[0];
rz(-1.7779508) q[1];
sx q[1];
rz(-1.9121871) q[1];
sx q[1];
rz(-1.8180064) q[1];
rz(-0.017756391) q[2];
sx q[2];
rz(-0.50928558) q[2];
sx q[2];
rz(-1.7094564) q[2];
rz(-1.2033403) q[3];
sx q[3];
rz(-2.8946579) q[3];
sx q[3];
rz(-2.2263087) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
