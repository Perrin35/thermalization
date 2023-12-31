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
rz(-2.354061) q[1];
sx q[1];
rz(-2.1252316) q[1];
sx q[1];
rz(-2.7273942) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83672268) q[0];
sx q[0];
rz(-2.6430336) q[0];
sx q[0];
rz(-2.3303633) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.094751058) q[2];
sx q[2];
rz(-2.13846) q[2];
sx q[2];
rz(-1.3324347) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.84250433) q[1];
sx q[1];
rz(-0.43049225) q[1];
sx q[1];
rz(1.4166142) q[1];
rz(-pi) q[2];
rz(0.26478404) q[3];
sx q[3];
rz(-1.4720819) q[3];
sx q[3];
rz(-2.7717154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9238613) q[2];
sx q[2];
rz(-1.2552746) q[2];
sx q[2];
rz(-0.031575354) q[2];
rz(1.8850373) q[3];
sx q[3];
rz(-2.6387408) q[3];
sx q[3];
rz(-0.52662915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3027705) q[0];
sx q[0];
rz(-1.4571723) q[0];
sx q[0];
rz(-0.1698499) q[0];
rz(-2.4376712) q[1];
sx q[1];
rz(-1.0715276) q[1];
sx q[1];
rz(0.53952113) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8115494) q[0];
sx q[0];
rz(-2.0199611) q[0];
sx q[0];
rz(-0.051785843) q[0];
x q[1];
rz(-0.72210724) q[2];
sx q[2];
rz(-1.2566084) q[2];
sx q[2];
rz(-2.9177641) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0454228) q[1];
sx q[1];
rz(-1.759937) q[1];
sx q[1];
rz(-1.063709) q[1];
rz(-pi) q[2];
rz(2.2803454) q[3];
sx q[3];
rz(-1.1303139) q[3];
sx q[3];
rz(1.4351821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.66723055) q[2];
sx q[2];
rz(-1.2381866) q[2];
sx q[2];
rz(0.24307069) q[2];
rz(2.4754751) q[3];
sx q[3];
rz(-0.56454286) q[3];
sx q[3];
rz(-1.8977785) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3617525) q[0];
sx q[0];
rz(-0.11479522) q[0];
sx q[0];
rz(0.4483805) q[0];
rz(-1.386863) q[1];
sx q[1];
rz(-1.153839) q[1];
sx q[1];
rz(-2.8853436) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2494333) q[0];
sx q[0];
rz(-1.8349577) q[0];
sx q[0];
rz(0.63412068) q[0];
rz(-2.3149812) q[2];
sx q[2];
rz(-0.6140784) q[2];
sx q[2];
rz(-0.6073063) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.2910034) q[1];
sx q[1];
rz(-2.2271419) q[1];
sx q[1];
rz(-2.7235051) q[1];
rz(-pi) q[2];
rz(-2.933421) q[3];
sx q[3];
rz(-2.9609207) q[3];
sx q[3];
rz(2.4908623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3391352) q[2];
sx q[2];
rz(-2.4352303) q[2];
sx q[2];
rz(-2.5668872) q[2];
rz(1.7859219) q[3];
sx q[3];
rz(-1.971743) q[3];
sx q[3];
rz(1.2333966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.213585) q[0];
sx q[0];
rz(-1.428823) q[0];
sx q[0];
rz(2.8821049) q[0];
rz(1.150594) q[1];
sx q[1];
rz(-1.3535627) q[1];
sx q[1];
rz(-2.4096699) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5808125) q[0];
sx q[0];
rz(-0.73045759) q[0];
sx q[0];
rz(-2.3985582) q[0];
rz(-3.0078997) q[2];
sx q[2];
rz(-2.4198654) q[2];
sx q[2];
rz(1.114811) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6432453) q[1];
sx q[1];
rz(-0.32532641) q[1];
sx q[1];
rz(0.64756067) q[1];
rz(-2.9210806) q[3];
sx q[3];
rz(-1.4482499) q[3];
sx q[3];
rz(1.8054655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6966454) q[2];
sx q[2];
rz(-1.8680633) q[2];
sx q[2];
rz(0.15110061) q[2];
rz(-0.54667306) q[3];
sx q[3];
rz(-2.0918545) q[3];
sx q[3];
rz(-2.7643519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5916409) q[0];
sx q[0];
rz(-1.0629835) q[0];
sx q[0];
rz(1.8792101) q[0];
rz(1.4683912) q[1];
sx q[1];
rz(-0.60931283) q[1];
sx q[1];
rz(-2.343822) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20194963) q[0];
sx q[0];
rz(-1.573283) q[0];
sx q[0];
rz(1.6909088) q[0];
x q[1];
rz(-0.30349489) q[2];
sx q[2];
rz(-1.7804838) q[2];
sx q[2];
rz(1.4022624) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3623912) q[1];
sx q[1];
rz(-0.63450846) q[1];
sx q[1];
rz(-1.4125376) q[1];
rz(-0.4425211) q[3];
sx q[3];
rz(-1.2708775) q[3];
sx q[3];
rz(0.72344852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.34565869) q[2];
sx q[2];
rz(-2.5107333) q[2];
sx q[2];
rz(0.30203715) q[2];
rz(1.9942412) q[3];
sx q[3];
rz(-1.6796422) q[3];
sx q[3];
rz(-2.5938477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(1.0844768) q[1];
sx q[1];
rz(-2.1613354) q[1];
sx q[1];
rz(-3.0715122) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4159211) q[0];
sx q[0];
rz(-1.7416746) q[0];
sx q[0];
rz(-2.9861949) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3577865) q[2];
sx q[2];
rz(-1.1597826) q[2];
sx q[2];
rz(2.59936) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.97092123) q[1];
sx q[1];
rz(-2.1340003) q[1];
sx q[1];
rz(1.3794273) q[1];
rz(1.4962247) q[3];
sx q[3];
rz(-1.6061022) q[3];
sx q[3];
rz(-1.9880184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8391116) q[2];
sx q[2];
rz(-1.182686) q[2];
sx q[2];
rz(0.62136674) q[2];
rz(-1.7012043) q[3];
sx q[3];
rz(-2.6337603) q[3];
sx q[3];
rz(-2.5682209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0550585) q[0];
sx q[0];
rz(-1.7250412) q[0];
sx q[0];
rz(0.85987464) q[0];
rz(1.2043918) q[1];
sx q[1];
rz(-0.87221968) q[1];
sx q[1];
rz(0.0079356114) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3646274) q[0];
sx q[0];
rz(-2.1205175) q[0];
sx q[0];
rz(-1.8718029) q[0];
rz(0.87848778) q[2];
sx q[2];
rz(-1.2307067) q[2];
sx q[2];
rz(-2.0932587) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.36564) q[1];
sx q[1];
rz(-2.0740293) q[1];
sx q[1];
rz(0.788049) q[1];
x q[2];
rz(-0.90585917) q[3];
sx q[3];
rz(-0.48478904) q[3];
sx q[3];
rz(-2.0743899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.69592151) q[2];
sx q[2];
rz(-1.9182703) q[2];
sx q[2];
rz(-0.41637862) q[2];
rz(-1.3683866) q[3];
sx q[3];
rz(-1.2973283) q[3];
sx q[3];
rz(-2.2369475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1798582) q[0];
sx q[0];
rz(-2.8680153) q[0];
sx q[0];
rz(0.36488786) q[0];
rz(-0.94003135) q[1];
sx q[1];
rz(-0.53661984) q[1];
sx q[1];
rz(1.6392802) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2973328) q[0];
sx q[0];
rz(-1.5607921) q[0];
sx q[0];
rz(-0.25015932) q[0];
rz(-pi) q[1];
rz(-1.0566696) q[2];
sx q[2];
rz(-2.3377315) q[2];
sx q[2];
rz(2.4664997) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8870526) q[1];
sx q[1];
rz(-0.92003838) q[1];
sx q[1];
rz(0.73144967) q[1];
rz(-pi) q[2];
rz(-3.0495166) q[3];
sx q[3];
rz(-1.7756988) q[3];
sx q[3];
rz(2.3931707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.49729785) q[2];
sx q[2];
rz(-0.50540322) q[2];
sx q[2];
rz(1.0650744) q[2];
rz(-2.8403357) q[3];
sx q[3];
rz(-1.4091636) q[3];
sx q[3];
rz(1.6206954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.168468) q[0];
sx q[0];
rz(-0.081806101) q[0];
sx q[0];
rz(0.43564963) q[0];
rz(1.7565953) q[1];
sx q[1];
rz(-0.47043097) q[1];
sx q[1];
rz(2.7246144) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0458841) q[0];
sx q[0];
rz(-2.4864712) q[0];
sx q[0];
rz(1.5089633) q[0];
rz(-pi) q[1];
rz(-2.4776393) q[2];
sx q[2];
rz(-1.10154) q[2];
sx q[2];
rz(-2.9479153) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3560564) q[1];
sx q[1];
rz(-2.0691263) q[1];
sx q[1];
rz(-2.9836125) q[1];
rz(-3.0005089) q[3];
sx q[3];
rz(-1.4420296) q[3];
sx q[3];
rz(-1.5878549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.10432648) q[2];
sx q[2];
rz(-1.6486847) q[2];
sx q[2];
rz(0.22582516) q[2];
rz(-0.2078235) q[3];
sx q[3];
rz(-0.72312975) q[3];
sx q[3];
rz(2.5549755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41480961) q[0];
sx q[0];
rz(-2.8746958) q[0];
sx q[0];
rz(-1.6171932) q[0];
rz(2.1879451) q[1];
sx q[1];
rz(-1.8667659) q[1];
sx q[1];
rz(-1.3226002) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1117301) q[0];
sx q[0];
rz(-2.7435281) q[0];
sx q[0];
rz(1.8882621) q[0];
rz(1.1216713) q[2];
sx q[2];
rz(-1.7977062) q[2];
sx q[2];
rz(0.34044468) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0304071) q[1];
sx q[1];
rz(-0.92799312) q[1];
sx q[1];
rz(-2.2305957) q[1];
rz(-0.91536509) q[3];
sx q[3];
rz(-1.4487106) q[3];
sx q[3];
rz(-1.9742427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7913197) q[2];
sx q[2];
rz(-2.0952756) q[2];
sx q[2];
rz(-1.0160758) q[2];
rz(1.2223876) q[3];
sx q[3];
rz(-0.17761579) q[3];
sx q[3];
rz(2.5861752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.9983457) q[0];
sx q[0];
rz(-2.2194942) q[0];
sx q[0];
rz(1.0794328) q[0];
rz(-1.3636419) q[1];
sx q[1];
rz(-1.2294055) q[1];
sx q[1];
rz(1.3235863) q[1];
rz(0.50921847) q[2];
sx q[2];
rz(-1.5794532) q[2];
sx q[2];
rz(3.0184359) q[2];
rz(0.090311269) q[3];
sx q[3];
rz(-1.3406546) q[3];
sx q[3];
rz(-1.8484074) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
