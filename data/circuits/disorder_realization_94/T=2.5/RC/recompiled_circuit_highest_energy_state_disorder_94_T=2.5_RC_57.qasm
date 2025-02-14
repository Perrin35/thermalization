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
rz(2.0512407) q[0];
sx q[0];
rz(-2.2418689) q[0];
sx q[0];
rz(-2.0119014) q[0];
rz(0.55454412) q[1];
sx q[1];
rz(-0.50538844) q[1];
sx q[1];
rz(-0.55355054) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3239297) q[0];
sx q[0];
rz(-0.50369064) q[0];
sx q[0];
rz(1.5445185) q[0];
x q[1];
rz(1.9808018) q[2];
sx q[2];
rz(-1.5799109) q[2];
sx q[2];
rz(1.3076919) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3379711) q[1];
sx q[1];
rz(-2.3052678) q[1];
sx q[1];
rz(2.5938889) q[1];
x q[2];
rz(-1.3517836) q[3];
sx q[3];
rz(-1.5288282) q[3];
sx q[3];
rz(-1.7445843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0022137) q[2];
sx q[2];
rz(-1.4993818) q[2];
sx q[2];
rz(1.1688983) q[2];
rz(-0.72541952) q[3];
sx q[3];
rz(-1.231671) q[3];
sx q[3];
rz(1.5861082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2296974) q[0];
sx q[0];
rz(-2.732369) q[0];
sx q[0];
rz(2.0912066) q[0];
rz(1.7627675) q[1];
sx q[1];
rz(-1.4070815) q[1];
sx q[1];
rz(-2.1398267) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9501099) q[0];
sx q[0];
rz(-1.313173) q[0];
sx q[0];
rz(0.10665032) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.839532) q[2];
sx q[2];
rz(-2.0141057) q[2];
sx q[2];
rz(-1.4745148) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6607501) q[1];
sx q[1];
rz(-2.6482812) q[1];
sx q[1];
rz(0.50220614) q[1];
rz(-pi) q[2];
rz(1.2779986) q[3];
sx q[3];
rz(-1.9778041) q[3];
sx q[3];
rz(-1.7990822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.072814552) q[2];
sx q[2];
rz(-0.38663703) q[2];
sx q[2];
rz(-1.6437423) q[2];
rz(2.1554598) q[3];
sx q[3];
rz(-2.17406) q[3];
sx q[3];
rz(-2.6223474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7943952) q[0];
sx q[0];
rz(-0.47413844) q[0];
sx q[0];
rz(-0.073632181) q[0];
rz(0.73078784) q[1];
sx q[1];
rz(-1.7981139) q[1];
sx q[1];
rz(-0.98339072) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7950566) q[0];
sx q[0];
rz(-1.7253526) q[0];
sx q[0];
rz(1.8554423) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.86536256) q[2];
sx q[2];
rz(-1.1956805) q[2];
sx q[2];
rz(0.6626216) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9485452) q[1];
sx q[1];
rz(-1.8406423) q[1];
sx q[1];
rz(0.13040925) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4296823) q[3];
sx q[3];
rz(-2.634438) q[3];
sx q[3];
rz(0.53559723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.694515) q[2];
sx q[2];
rz(-2.647001) q[2];
sx q[2];
rz(3.0818603) q[2];
rz(-0.51074243) q[3];
sx q[3];
rz(-2.7849438) q[3];
sx q[3];
rz(1.3753148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.315595) q[0];
sx q[0];
rz(-0.70846486) q[0];
sx q[0];
rz(0.84061709) q[0];
rz(0.21545848) q[1];
sx q[1];
rz(-1.0865728) q[1];
sx q[1];
rz(-0.1159018) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.602086) q[0];
sx q[0];
rz(-1.2316252) q[0];
sx q[0];
rz(2.8643292) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8491146) q[2];
sx q[2];
rz(-1.3693606) q[2];
sx q[2];
rz(-0.62244994) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.85899587) q[1];
sx q[1];
rz(-2.1934185) q[1];
sx q[1];
rz(-2.7864561) q[1];
rz(0.41036116) q[3];
sx q[3];
rz(-1.6660353) q[3];
sx q[3];
rz(1.8382753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0243715) q[2];
sx q[2];
rz(-1.2933967) q[2];
sx q[2];
rz(0.8563861) q[2];
rz(0.74445009) q[3];
sx q[3];
rz(-0.88862935) q[3];
sx q[3];
rz(0.3652679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.982548) q[0];
sx q[0];
rz(-0.58369517) q[0];
sx q[0];
rz(2.1481376) q[0];
rz(-1.1194725) q[1];
sx q[1];
rz(-1.9529724) q[1];
sx q[1];
rz(-3.1192034) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36528811) q[0];
sx q[0];
rz(-0.1488221) q[0];
sx q[0];
rz(-2.5689199) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55565067) q[2];
sx q[2];
rz(-2.0731063) q[2];
sx q[2];
rz(-2.7567418) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7931079) q[1];
sx q[1];
rz(-2.1505618) q[1];
sx q[1];
rz(2.0391885) q[1];
rz(-pi) q[2];
rz(2.0517572) q[3];
sx q[3];
rz(-2.2053455) q[3];
sx q[3];
rz(-0.1635199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.28919724) q[2];
sx q[2];
rz(-0.22326938) q[2];
sx q[2];
rz(1.7247464) q[2];
rz(2.1220186) q[3];
sx q[3];
rz(-2.151078) q[3];
sx q[3];
rz(0.31653252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9483017) q[0];
sx q[0];
rz(-3.1315849) q[0];
sx q[0];
rz(-0.77578068) q[0];
rz(-1.4099482) q[1];
sx q[1];
rz(-1.3103139) q[1];
sx q[1];
rz(-1.5595248) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0163859) q[0];
sx q[0];
rz(-1.6324955) q[0];
sx q[0];
rz(0.030304219) q[0];
rz(2.1072793) q[2];
sx q[2];
rz(-0.63648495) q[2];
sx q[2];
rz(0.45490593) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0770862) q[1];
sx q[1];
rz(-0.48708579) q[1];
sx q[1];
rz(1.8555831) q[1];
x q[2];
rz(-2.6234145) q[3];
sx q[3];
rz(-1.2478531) q[3];
sx q[3];
rz(-0.87080982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5887624) q[2];
sx q[2];
rz(-1.4344119) q[2];
sx q[2];
rz(0.32858783) q[2];
rz(-1.5821247) q[3];
sx q[3];
rz(-2.4170473) q[3];
sx q[3];
rz(-2.4174378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5174823) q[0];
sx q[0];
rz(-0.8172074) q[0];
sx q[0];
rz(2.4716614) q[0];
rz(-1.6239369) q[1];
sx q[1];
rz(-1.0847963) q[1];
sx q[1];
rz(-2.0673015) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6321323) q[0];
sx q[0];
rz(-1.7863417) q[0];
sx q[0];
rz(-2.4556405) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7781717) q[2];
sx q[2];
rz(-1.3989903) q[2];
sx q[2];
rz(-1.2762783) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.049557471) q[1];
sx q[1];
rz(-1.2566031) q[1];
sx q[1];
rz(-0.097196984) q[1];
rz(-pi) q[2];
rz(-0.89364918) q[3];
sx q[3];
rz(-1.4010249) q[3];
sx q[3];
rz(2.5285965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.39277789) q[2];
sx q[2];
rz(-2.1608976) q[2];
sx q[2];
rz(-1.7559715) q[2];
rz(-0.70220002) q[3];
sx q[3];
rz(-2.8148013) q[3];
sx q[3];
rz(-1.1087815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4362519) q[0];
sx q[0];
rz(-1.5441283) q[0];
sx q[0];
rz(2.3633549) q[0];
rz(-1.8477919) q[1];
sx q[1];
rz(-1.8840645) q[1];
sx q[1];
rz(-0.055880849) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.125408) q[0];
sx q[0];
rz(-1.2795537) q[0];
sx q[0];
rz(-1.6470558) q[0];
rz(2.4498522) q[2];
sx q[2];
rz(-2.8517339) q[2];
sx q[2];
rz(-0.65164072) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3371007) q[1];
sx q[1];
rz(-0.40834112) q[1];
sx q[1];
rz(-0.33894811) q[1];
rz(-pi) q[2];
rz(-2.3010766) q[3];
sx q[3];
rz(-1.9728246) q[3];
sx q[3];
rz(-0.4268643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.087626584) q[2];
sx q[2];
rz(-2.1199333) q[2];
sx q[2];
rz(2.450024) q[2];
rz(-0.21928731) q[3];
sx q[3];
rz(-1.4838452) q[3];
sx q[3];
rz(-1.1933491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9241065) q[0];
sx q[0];
rz(-3.1052488) q[0];
sx q[0];
rz(-2.0935667) q[0];
rz(-2.3955087) q[1];
sx q[1];
rz(-1.8522976) q[1];
sx q[1];
rz(-2.7229436) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9384883) q[0];
sx q[0];
rz(-3.1230429) q[0];
sx q[0];
rz(2.0880329) q[0];
x q[1];
rz(0.92709686) q[2];
sx q[2];
rz(-1.3852556) q[2];
sx q[2];
rz(-0.58840051) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9050952) q[1];
sx q[1];
rz(-1.7069382) q[1];
sx q[1];
rz(-0.13827583) q[1];
rz(-3.0935086) q[3];
sx q[3];
rz(-2.2732031) q[3];
sx q[3];
rz(2.4772205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5082385) q[2];
sx q[2];
rz(-0.60428047) q[2];
sx q[2];
rz(-2.5227127) q[2];
rz(-0.45281705) q[3];
sx q[3];
rz(-1.6524977) q[3];
sx q[3];
rz(0.49310163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.598269) q[0];
sx q[0];
rz(-3.0040574) q[0];
sx q[0];
rz(1.2602873) q[0];
rz(-2.8434985) q[1];
sx q[1];
rz(-2.1838078) q[1];
sx q[1];
rz(1.4990998) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74841046) q[0];
sx q[0];
rz(-1.5934738) q[0];
sx q[0];
rz(-2.3298954) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7146729) q[2];
sx q[2];
rz(-1.292079) q[2];
sx q[2];
rz(2.4867694) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.11718845) q[1];
sx q[1];
rz(-2.9147101) q[1];
sx q[1];
rz(1.4972761) q[1];
rz(-pi) q[2];
rz(-1.9630249) q[3];
sx q[3];
rz(-1.4334599) q[3];
sx q[3];
rz(1.3016421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5160211) q[2];
sx q[2];
rz(-1.2837774) q[2];
sx q[2];
rz(2.9594291) q[2];
rz(-0.22570172) q[3];
sx q[3];
rz(-1.00939) q[3];
sx q[3];
rz(-1.3753447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7094649) q[0];
sx q[0];
rz(-1.3668677) q[0];
sx q[0];
rz(0.9912542) q[0];
rz(1.2200914) q[1];
sx q[1];
rz(-2.5972001) q[1];
sx q[1];
rz(-0.77144815) q[1];
rz(1.3924072) q[2];
sx q[2];
rz(-1.6606944) q[2];
sx q[2];
rz(0.14002945) q[2];
rz(2.3301043) q[3];
sx q[3];
rz(-0.93899721) q[3];
sx q[3];
rz(-0.39427857) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
