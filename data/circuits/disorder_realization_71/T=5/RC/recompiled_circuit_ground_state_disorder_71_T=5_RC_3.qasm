OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7856287) q[0];
sx q[0];
rz(-0.2956737) q[0];
sx q[0];
rz(0.41391882) q[0];
rz(3.4105372) q[1];
sx q[1];
rz(6.6528448) q[1];
sx q[1];
rz(11.289968) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7463425) q[0];
sx q[0];
rz(-1.479404) q[0];
sx q[0];
rz(1.6618098) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4914614) q[2];
sx q[2];
rz(-2.4717307) q[2];
sx q[2];
rz(1.2487433) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2416083) q[1];
sx q[1];
rz(-1.2869724) q[1];
sx q[1];
rz(0.79197661) q[1];
rz(-pi) q[2];
rz(-2.5475575) q[3];
sx q[3];
rz(-1.4184991) q[3];
sx q[3];
rz(1.9826979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.74497491) q[2];
sx q[2];
rz(-1.9048385) q[2];
sx q[2];
rz(1.2006987) q[2];
rz(-0.68136224) q[3];
sx q[3];
rz(-1.5228289) q[3];
sx q[3];
rz(2.7549226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6354527) q[0];
sx q[0];
rz(-2.8400087) q[0];
sx q[0];
rz(2.8801081) q[0];
rz(1.6188949) q[1];
sx q[1];
rz(-2.6324184) q[1];
sx q[1];
rz(1.8656628) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0715213) q[0];
sx q[0];
rz(-2.1198648) q[0];
sx q[0];
rz(1.4577876) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66565565) q[2];
sx q[2];
rz(-2.3710069) q[2];
sx q[2];
rz(1.6536755) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9918242) q[1];
sx q[1];
rz(-2.5006378) q[1];
sx q[1];
rz(0.561552) q[1];
rz(3.1072282) q[3];
sx q[3];
rz(-1.7196764) q[3];
sx q[3];
rz(1.9158165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2168938) q[2];
sx q[2];
rz(-2.108722) q[2];
sx q[2];
rz(-0.72803289) q[2];
rz(-1.3701471) q[3];
sx q[3];
rz(-0.61714211) q[3];
sx q[3];
rz(-0.035577687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.034123357) q[0];
sx q[0];
rz(-2.0296445) q[0];
sx q[0];
rz(2.4566101) q[0];
rz(-2.0383539) q[1];
sx q[1];
rz(-2.5220242) q[1];
sx q[1];
rz(1.0122976) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023803614) q[0];
sx q[0];
rz(-0.53142457) q[0];
sx q[0];
rz(0.84367911) q[0];
rz(2.9673789) q[2];
sx q[2];
rz(-1.3784587) q[2];
sx q[2];
rz(1.9680374) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5033414) q[1];
sx q[1];
rz(-1.4443399) q[1];
sx q[1];
rz(-2.9328037) q[1];
x q[2];
rz(-0.52148444) q[3];
sx q[3];
rz(-2.3064724) q[3];
sx q[3];
rz(0.84671742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.181902) q[2];
sx q[2];
rz(-1.5896553) q[2];
sx q[2];
rz(2.336179) q[2];
rz(2.3253333) q[3];
sx q[3];
rz(-2.5160242) q[3];
sx q[3];
rz(-0.094154112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1212946) q[0];
sx q[0];
rz(-2.2327406) q[0];
sx q[0];
rz(-0.62776172) q[0];
rz(-0.10920564) q[1];
sx q[1];
rz(-0.86995482) q[1];
sx q[1];
rz(-2.3039718) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35204404) q[0];
sx q[0];
rz(-1.0800011) q[0];
sx q[0];
rz(1.646203) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2455279) q[2];
sx q[2];
rz(-0.62042371) q[2];
sx q[2];
rz(-1.9577946) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.301103) q[1];
sx q[1];
rz(-0.87574848) q[1];
sx q[1];
rz(-1.3962505) q[1];
rz(-1.8786976) q[3];
sx q[3];
rz(-1.2114204) q[3];
sx q[3];
rz(-1.9482159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3179021) q[2];
sx q[2];
rz(-0.48419848) q[2];
sx q[2];
rz(2.79706) q[2];
rz(-1.7349998) q[3];
sx q[3];
rz(-2.78648) q[3];
sx q[3];
rz(0.047903456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59871167) q[0];
sx q[0];
rz(-2.5051835) q[0];
sx q[0];
rz(2.6357546) q[0];
rz(2.089031) q[1];
sx q[1];
rz(-2.274175) q[1];
sx q[1];
rz(0.94598407) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.036412843) q[0];
sx q[0];
rz(-2.7650021) q[0];
sx q[0];
rz(0.74684493) q[0];
rz(-0.35752488) q[2];
sx q[2];
rz(-1.4357743) q[2];
sx q[2];
rz(2.6808878) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.16984384) q[1];
sx q[1];
rz(-1.1976114) q[1];
sx q[1];
rz(-0.62974522) q[1];
x q[2];
rz(1.9993773) q[3];
sx q[3];
rz(-0.16332291) q[3];
sx q[3];
rz(-1.6960916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0537009) q[2];
sx q[2];
rz(-1.0975857) q[2];
sx q[2];
rz(1.5637448) q[2];
rz(2.0569233) q[3];
sx q[3];
rz(-2.2659149) q[3];
sx q[3];
rz(-0.55115551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56187335) q[0];
sx q[0];
rz(-2.5358574) q[0];
sx q[0];
rz(0.24000034) q[0];
rz(-0.9217681) q[1];
sx q[1];
rz(-1.0488989) q[1];
sx q[1];
rz(1.1704495) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76239785) q[0];
sx q[0];
rz(-2.1403011) q[0];
sx q[0];
rz(1.5194917) q[0];
rz(-1.6977152) q[2];
sx q[2];
rz(-1.2170212) q[2];
sx q[2];
rz(-0.60777174) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.45615087) q[1];
sx q[1];
rz(-0.69760453) q[1];
sx q[1];
rz(-0.82932034) q[1];
rz(-pi) q[2];
rz(-1.7975054) q[3];
sx q[3];
rz(-2.7313958) q[3];
sx q[3];
rz(2.6997379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.74448284) q[2];
sx q[2];
rz(-2.3220671) q[2];
sx q[2];
rz(-0.58103713) q[2];
rz(2.2033612) q[3];
sx q[3];
rz(-0.56648985) q[3];
sx q[3];
rz(-0.70403045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0051603) q[0];
sx q[0];
rz(-0.67661023) q[0];
sx q[0];
rz(-0.50049472) q[0];
rz(-2.9035134) q[1];
sx q[1];
rz(-1.6903189) q[1];
sx q[1];
rz(-2.4949825) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71599267) q[0];
sx q[0];
rz(-1.5690205) q[0];
sx q[0];
rz(0.034816936) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2836566) q[2];
sx q[2];
rz(-0.11936766) q[2];
sx q[2];
rz(-2.067208) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2518066) q[1];
sx q[1];
rz(-1.4507035) q[1];
sx q[1];
rz(-2.0577862) q[1];
rz(-pi) q[2];
rz(1.7342689) q[3];
sx q[3];
rz(-0.36296926) q[3];
sx q[3];
rz(-0.27856871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.87767345) q[2];
sx q[2];
rz(-0.71920243) q[2];
sx q[2];
rz(0.70362299) q[2];
rz(-1.2879114) q[3];
sx q[3];
rz(-1.0602919) q[3];
sx q[3];
rz(-0.50584403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8877761) q[0];
sx q[0];
rz(-1.1340589) q[0];
sx q[0];
rz(1.8164841) q[0];
rz(-0.77249402) q[1];
sx q[1];
rz(-1.1224116) q[1];
sx q[1];
rz(-0.17556369) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2529375) q[0];
sx q[0];
rz(-2.3217259) q[0];
sx q[0];
rz(1.3142513) q[0];
rz(2.2261691) q[2];
sx q[2];
rz(-0.88214126) q[2];
sx q[2];
rz(1.6885533) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.23607924) q[1];
sx q[1];
rz(-1.7100952) q[1];
sx q[1];
rz(1.1054611) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14388559) q[3];
sx q[3];
rz(-0.7215313) q[3];
sx q[3];
rz(0.084621457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0385711) q[2];
sx q[2];
rz(-2.0671637) q[2];
sx q[2];
rz(-1.0938905) q[2];
rz(-1.966018) q[3];
sx q[3];
rz(-2.3837377) q[3];
sx q[3];
rz(0.28483835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4232101) q[0];
sx q[0];
rz(-0.045561401) q[0];
sx q[0];
rz(1.3257931) q[0];
rz(2.6153053) q[1];
sx q[1];
rz(-0.86295366) q[1];
sx q[1];
rz(2.3525499) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5518582) q[0];
sx q[0];
rz(-1.6673256) q[0];
sx q[0];
rz(-1.4087229) q[0];
x q[1];
rz(0.61179917) q[2];
sx q[2];
rz(-0.42610301) q[2];
sx q[2];
rz(-0.95624051) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.583786) q[1];
sx q[1];
rz(-0.72804385) q[1];
sx q[1];
rz(-1.0398204) q[1];
rz(2.3655543) q[3];
sx q[3];
rz(-2.0774842) q[3];
sx q[3];
rz(2.6069178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9290756) q[2];
sx q[2];
rz(-0.93364659) q[2];
sx q[2];
rz(0.90544256) q[2];
rz(2.2255911) q[3];
sx q[3];
rz(-1.6986877) q[3];
sx q[3];
rz(1.0441095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4395897) q[0];
sx q[0];
rz(-2.3083394) q[0];
sx q[0];
rz(0.92779094) q[0];
rz(-2.0480428) q[1];
sx q[1];
rz(-0.55892006) q[1];
sx q[1];
rz(-3.0082026) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0722712) q[0];
sx q[0];
rz(-1.7786553) q[0];
sx q[0];
rz(0.29671191) q[0];
x q[1];
rz(0.20602823) q[2];
sx q[2];
rz(-1.0090172) q[2];
sx q[2];
rz(-2.5453407) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9815753) q[1];
sx q[1];
rz(-0.80343738) q[1];
sx q[1];
rz(2.395242) q[1];
rz(-pi) q[2];
rz(1.1130775) q[3];
sx q[3];
rz(-1.7594171) q[3];
sx q[3];
rz(-0.45841226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6371969) q[2];
sx q[2];
rz(-1.7375676) q[2];
sx q[2];
rz(-1.8782328) q[2];
rz(1.0697621) q[3];
sx q[3];
rz(-2.1061149) q[3];
sx q[3];
rz(-0.076816471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99061154) q[0];
sx q[0];
rz(-2.4438416) q[0];
sx q[0];
rz(-1.7846815) q[0];
rz(-1.6928584) q[1];
sx q[1];
rz(-2.3285463) q[1];
sx q[1];
rz(2.9978233) q[1];
rz(-2.095456) q[2];
sx q[2];
rz(-1.2992278) q[2];
sx q[2];
rz(-2.3067731) q[2];
rz(-1.2273905) q[3];
sx q[3];
rz(-1.3385942) q[3];
sx q[3];
rz(-0.75029324) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
