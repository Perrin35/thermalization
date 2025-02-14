OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.4177999) q[0];
sx q[0];
rz(-1.1136709) q[0];
sx q[0];
rz(-0.96473515) q[0];
rz(-0.073702987) q[1];
sx q[1];
rz(-0.6414203) q[1];
sx q[1];
rz(-1.4944271) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9091563) q[0];
sx q[0];
rz(-2.3415544) q[0];
sx q[0];
rz(2.590488) q[0];
rz(-pi) q[1];
rz(3.037248) q[2];
sx q[2];
rz(-1.8868703) q[2];
sx q[2];
rz(0.8546654) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9508178) q[1];
sx q[1];
rz(-0.72682805) q[1];
sx q[1];
rz(-0.72587691) q[1];
rz(-pi) q[2];
rz(-2.1001746) q[3];
sx q[3];
rz(-1.6106859) q[3];
sx q[3];
rz(-0.94464579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.73704314) q[2];
sx q[2];
rz(-1.6208384) q[2];
sx q[2];
rz(2.9250195) q[2];
rz(2.6259322) q[3];
sx q[3];
rz(-2.0629081) q[3];
sx q[3];
rz(-3.0576341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11230042) q[0];
sx q[0];
rz(-3.0700505) q[0];
sx q[0];
rz(1.4805502) q[0];
rz(2.5484565) q[1];
sx q[1];
rz(-2.1430404) q[1];
sx q[1];
rz(0.28233972) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7337228) q[0];
sx q[0];
rz(-1.5998915) q[0];
sx q[0];
rz(0.12761527) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7690053) q[2];
sx q[2];
rz(-1.7019801) q[2];
sx q[2];
rz(-0.030156915) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8806539) q[1];
sx q[1];
rz(-1.0323413) q[1];
sx q[1];
rz(1.4101342) q[1];
x q[2];
rz(2.2766791) q[3];
sx q[3];
rz(-2.5967715) q[3];
sx q[3];
rz(-0.13710216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4234408) q[2];
sx q[2];
rz(-0.56090063) q[2];
sx q[2];
rz(2.1924428) q[2];
rz(-0.078350457) q[3];
sx q[3];
rz(-1.6716985) q[3];
sx q[3];
rz(1.9115492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3345728) q[0];
sx q[0];
rz(-3.0834575) q[0];
sx q[0];
rz(-0.19244254) q[0];
rz(-0.9737393) q[1];
sx q[1];
rz(-2.1111919) q[1];
sx q[1];
rz(1.8691501) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21308831) q[0];
sx q[0];
rz(-2.4165389) q[0];
sx q[0];
rz(-0.37037767) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3272304) q[2];
sx q[2];
rz(-2.3911016) q[2];
sx q[2];
rz(-1.0433973) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7682926) q[1];
sx q[1];
rz(-1.8197316) q[1];
sx q[1];
rz(1.8087923) q[1];
rz(-1.1149241) q[3];
sx q[3];
rz(-1.6041363) q[3];
sx q[3];
rz(0.10297266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.56866139) q[2];
sx q[2];
rz(-0.91721407) q[2];
sx q[2];
rz(0.60079637) q[2];
rz(2.3066547) q[3];
sx q[3];
rz(-2.2227414) q[3];
sx q[3];
rz(0.44343534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7610953) q[0];
sx q[0];
rz(-2.0090071) q[0];
sx q[0];
rz(-1.4196716) q[0];
rz(0.28600606) q[1];
sx q[1];
rz(-2.0068469) q[1];
sx q[1];
rz(2.6223415) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38558593) q[0];
sx q[0];
rz(-2.950188) q[0];
sx q[0];
rz(0.49218224) q[0];
x q[1];
rz(-0.51651603) q[2];
sx q[2];
rz(-2.9331547) q[2];
sx q[2];
rz(-1.0112273) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9870465) q[1];
sx q[1];
rz(-2.0813165) q[1];
sx q[1];
rz(-1.4922008) q[1];
rz(1.5386228) q[3];
sx q[3];
rz(-1.0250499) q[3];
sx q[3];
rz(-1.7189027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9754192) q[2];
sx q[2];
rz(-2.4611046) q[2];
sx q[2];
rz(2.5892881) q[2];
rz(2.35899) q[3];
sx q[3];
rz(-1.3026594) q[3];
sx q[3];
rz(2.9569614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(0.85823378) q[0];
sx q[0];
rz(-2.5193546) q[0];
sx q[0];
rz(-1.8587814) q[0];
rz(1.917631) q[1];
sx q[1];
rz(-1.6061648) q[1];
sx q[1];
rz(0.70972365) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2252475) q[0];
sx q[0];
rz(-3.1053251) q[0];
sx q[0];
rz(-0.64860089) q[0];
x q[1];
rz(3.100509) q[2];
sx q[2];
rz(-2.6729995) q[2];
sx q[2];
rz(-3.1236609) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.57576631) q[1];
sx q[1];
rz(-1.499956) q[1];
sx q[1];
rz(1.2567568) q[1];
rz(-pi) q[2];
rz(1.1790347) q[3];
sx q[3];
rz(-0.44322792) q[3];
sx q[3];
rz(2.0947651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.84511406) q[2];
sx q[2];
rz(-0.60135403) q[2];
sx q[2];
rz(-0.55310407) q[2];
rz(2.3002355) q[3];
sx q[3];
rz(-2.0935121) q[3];
sx q[3];
rz(-2.0060189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3784921) q[0];
sx q[0];
rz(-1.0569514) q[0];
sx q[0];
rz(0.77504778) q[0];
rz(-1.4729602) q[1];
sx q[1];
rz(-0.62455636) q[1];
sx q[1];
rz(-0.83736173) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.178029) q[0];
sx q[0];
rz(-0.63321165) q[0];
sx q[0];
rz(-1.4073611) q[0];
rz(-2.6858575) q[2];
sx q[2];
rz(-1.4066469) q[2];
sx q[2];
rz(-1.8388621) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.130116) q[1];
sx q[1];
rz(-1.1728047) q[1];
sx q[1];
rz(0.46020646) q[1];
rz(1.0488408) q[3];
sx q[3];
rz(-2.6673954) q[3];
sx q[3];
rz(-1.3124581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2396635) q[2];
sx q[2];
rz(-1.9749125) q[2];
sx q[2];
rz(0.76266328) q[2];
rz(-0.20721063) q[3];
sx q[3];
rz(-1.8800294) q[3];
sx q[3];
rz(0.60780779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36829456) q[0];
sx q[0];
rz(-0.93450707) q[0];
sx q[0];
rz(1.4542798) q[0];
rz(-1.2794718) q[1];
sx q[1];
rz(-0.8431294) q[1];
sx q[1];
rz(-1.4131193) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9521546) q[0];
sx q[0];
rz(-1.7725774) q[0];
sx q[0];
rz(-0.39315572) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83820512) q[2];
sx q[2];
rz(-1.1229441) q[2];
sx q[2];
rz(-3.107576) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3154732) q[1];
sx q[1];
rz(-0.54448381) q[1];
sx q[1];
rz(0.61750169) q[1];
rz(1.132148) q[3];
sx q[3];
rz(-1.8941714) q[3];
sx q[3];
rz(-1.1667337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1450682) q[2];
sx q[2];
rz(-2.2729496) q[2];
sx q[2];
rz(0.16150148) q[2];
rz(0.7343556) q[3];
sx q[3];
rz(-0.86638325) q[3];
sx q[3];
rz(1.3041147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9228009) q[0];
sx q[0];
rz(-2.6113593) q[0];
sx q[0];
rz(2.9042322) q[0];
rz(2.3573719) q[1];
sx q[1];
rz(-1.7796681) q[1];
sx q[1];
rz(-1.721419) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.500017) q[0];
sx q[0];
rz(-0.58860129) q[0];
sx q[0];
rz(-2.5831548) q[0];
rz(-1.3018441) q[2];
sx q[2];
rz(-0.70640817) q[2];
sx q[2];
rz(3.067467) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8964149) q[1];
sx q[1];
rz(-0.62508622) q[1];
sx q[1];
rz(0.51621373) q[1];
rz(-pi) q[2];
x q[2];
rz(0.5077614) q[3];
sx q[3];
rz(-1.1089033) q[3];
sx q[3];
rz(-2.2162102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.64804849) q[2];
sx q[2];
rz(-1.8914787) q[2];
sx q[2];
rz(1.9165215) q[2];
rz(-1.8152292) q[3];
sx q[3];
rz(-2.1882961) q[3];
sx q[3];
rz(1.1939322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61106435) q[0];
sx q[0];
rz(-3.0623797) q[0];
sx q[0];
rz(0.63842574) q[0];
rz(-0.05052677) q[1];
sx q[1];
rz(-1.2082929) q[1];
sx q[1];
rz(-2.5343177) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3359708) q[0];
sx q[0];
rz(-0.50862776) q[0];
sx q[0];
rz(-0.87368272) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4953682) q[2];
sx q[2];
rz(-0.68623073) q[2];
sx q[2];
rz(-2.8099912) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.0083976204) q[1];
sx q[1];
rz(-1.5004769) q[1];
sx q[1];
rz(-1.8969055) q[1];
x q[2];
rz(1.7631986) q[3];
sx q[3];
rz(-1.3891313) q[3];
sx q[3];
rz(2.3013129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.66024572) q[2];
sx q[2];
rz(-1.9504184) q[2];
sx q[2];
rz(2.353239) q[2];
rz(-0.83769074) q[3];
sx q[3];
rz(-0.82101429) q[3];
sx q[3];
rz(1.8590417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7343219) q[0];
sx q[0];
rz(-2.4096074) q[0];
sx q[0];
rz(0.5087854) q[0];
rz(-1.58163) q[1];
sx q[1];
rz(-2.1211233) q[1];
sx q[1];
rz(-0.4090974) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1412774) q[0];
sx q[0];
rz(-2.8594236) q[0];
sx q[0];
rz(-0.98833584) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81014388) q[2];
sx q[2];
rz(-0.64254566) q[2];
sx q[2];
rz(1.0660397) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.075114014) q[1];
sx q[1];
rz(-1.7914322) q[1];
sx q[1];
rz(1.5599987) q[1];
rz(0.76830058) q[3];
sx q[3];
rz(-0.51737201) q[3];
sx q[3];
rz(2.9381772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.540655) q[2];
sx q[2];
rz(-1.5254075) q[2];
sx q[2];
rz(1.2826922) q[2];
rz(3.1278074) q[3];
sx q[3];
rz(-1.1252334) q[3];
sx q[3];
rz(1.8085326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29671058) q[0];
sx q[0];
rz(-1.7322576) q[0];
sx q[0];
rz(-2.5175293) q[0];
rz(-0.86112549) q[1];
sx q[1];
rz(-2.594941) q[1];
sx q[1];
rz(-3.0088967) q[1];
rz(1.672198) q[2];
sx q[2];
rz(-1.8894926) q[2];
sx q[2];
rz(2.8429902) q[2];
rz(-0.80879296) q[3];
sx q[3];
rz(-0.8142796) q[3];
sx q[3];
rz(0.7249226) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
