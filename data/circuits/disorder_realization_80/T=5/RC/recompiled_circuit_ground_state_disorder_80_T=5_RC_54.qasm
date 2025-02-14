OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.6516946) q[0];
sx q[0];
rz(-0.88180056) q[0];
sx q[0];
rz(-0.14468004) q[0];
rz(-2.7235003) q[1];
sx q[1];
rz(-3.0378208) q[1];
sx q[1];
rz(1.6296847) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3016925) q[0];
sx q[0];
rz(-1.4524676) q[0];
sx q[0];
rz(-2.0685643) q[0];
rz(-0.73349632) q[2];
sx q[2];
rz(-1.6111177) q[2];
sx q[2];
rz(2.2625543) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.50070333) q[1];
sx q[1];
rz(-1.2391866) q[1];
sx q[1];
rz(1.9808116) q[1];
rz(-pi) q[2];
rz(-0.51201323) q[3];
sx q[3];
rz(-1.1992992) q[3];
sx q[3];
rz(-3.01513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8251553) q[2];
sx q[2];
rz(-1.2707767) q[2];
sx q[2];
rz(2.4862508) q[2];
rz(-1.3842899) q[3];
sx q[3];
rz(-0.96450788) q[3];
sx q[3];
rz(0.82096076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4493988) q[0];
sx q[0];
rz(-1.3757502) q[0];
sx q[0];
rz(-0.50814116) q[0];
rz(-1.3805768) q[1];
sx q[1];
rz(-0.5813798) q[1];
sx q[1];
rz(1.9320206) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0469477) q[0];
sx q[0];
rz(-0.98804615) q[0];
sx q[0];
rz(-2.9215982) q[0];
x q[1];
rz(1.5884933) q[2];
sx q[2];
rz(-1.624384) q[2];
sx q[2];
rz(0.45547661) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.74871127) q[1];
sx q[1];
rz(-1.4429174) q[1];
sx q[1];
rz(-2.676968) q[1];
rz(-pi) q[2];
rz(-2.1749635) q[3];
sx q[3];
rz(-1.2467017) q[3];
sx q[3];
rz(1.5965727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2035344) q[2];
sx q[2];
rz(-1.7455696) q[2];
sx q[2];
rz(-0.61689845) q[2];
rz(-1.5361702) q[3];
sx q[3];
rz(-1.0441531) q[3];
sx q[3];
rz(2.6507586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1598375) q[0];
sx q[0];
rz(-1.6438537) q[0];
sx q[0];
rz(-0.72506654) q[0];
rz(1.8720576) q[1];
sx q[1];
rz(-0.67765647) q[1];
sx q[1];
rz(2.1824172) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1364258) q[0];
sx q[0];
rz(-1.1014897) q[0];
sx q[0];
rz(0.43715663) q[0];
rz(-pi) q[1];
rz(-2.3940635) q[2];
sx q[2];
rz(-0.5420712) q[2];
sx q[2];
rz(-0.67491787) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7669285) q[1];
sx q[1];
rz(-0.21682021) q[1];
sx q[1];
rz(-1.0105074) q[1];
x q[2];
rz(1.3884344) q[3];
sx q[3];
rz(-1.7373457) q[3];
sx q[3];
rz(-2.9349851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3265257) q[2];
sx q[2];
rz(-1.2667789) q[2];
sx q[2];
rz(0.25700021) q[2];
rz(-1.0810931) q[3];
sx q[3];
rz(-0.5210146) q[3];
sx q[3];
rz(-1.5628373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0839888) q[0];
sx q[0];
rz(-0.27844089) q[0];
sx q[0];
rz(1.1778911) q[0];
rz(0.15011694) q[1];
sx q[1];
rz(-0.53214407) q[1];
sx q[1];
rz(1.5163126) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1105791) q[0];
sx q[0];
rz(-2.7000801) q[0];
sx q[0];
rz(2.2792321) q[0];
x q[1];
rz(2.7706835) q[2];
sx q[2];
rz(-2.264655) q[2];
sx q[2];
rz(2.3313076) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.0016370694) q[1];
sx q[1];
rz(-1.4712584) q[1];
sx q[1];
rz(0.32685728) q[1];
x q[2];
rz(-1.7606973) q[3];
sx q[3];
rz(-2.2058124) q[3];
sx q[3];
rz(0.39228016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2991422) q[2];
sx q[2];
rz(-2.5717042) q[2];
sx q[2];
rz(-2.2184856) q[2];
rz(-2.182377) q[3];
sx q[3];
rz(-2.1289289) q[3];
sx q[3];
rz(-2.925351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(0.32090309) q[0];
sx q[0];
rz(-2.237759) q[0];
sx q[0];
rz(0.88291105) q[0];
rz(1.2954905) q[1];
sx q[1];
rz(-0.77406445) q[1];
sx q[1];
rz(2.6466218) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2705602) q[0];
sx q[0];
rz(-0.43631662) q[0];
sx q[0];
rz(0.83026921) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60762824) q[2];
sx q[2];
rz(-2.3403185) q[2];
sx q[2];
rz(-2.5382588) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0985097) q[1];
sx q[1];
rz(-1.483023) q[1];
sx q[1];
rz(2.4721774) q[1];
rz(-pi) q[2];
rz(0.99397387) q[3];
sx q[3];
rz(-0.60189542) q[3];
sx q[3];
rz(0.26443297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3317269) q[2];
sx q[2];
rz(-0.11955424) q[2];
sx q[2];
rz(1.2681819) q[2];
rz(3.0590893) q[3];
sx q[3];
rz(-1.637633) q[3];
sx q[3];
rz(0.65703195) q[3];
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
rz(-1.4853915) q[0];
sx q[0];
rz(-2.9347561) q[0];
sx q[0];
rz(-1.50151) q[0];
rz(-2.0791176) q[1];
sx q[1];
rz(-1.1524009) q[1];
sx q[1];
rz(-1.9662439) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51477434) q[0];
sx q[0];
rz(-1.2668224) q[0];
sx q[0];
rz(1.8016812) q[0];
rz(-0.34890596) q[2];
sx q[2];
rz(-2.4173173) q[2];
sx q[2];
rz(1.6394212) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2986002) q[1];
sx q[1];
rz(-1.6879884) q[1];
sx q[1];
rz(0.25335626) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3812416) q[3];
sx q[3];
rz(-1.3521128) q[3];
sx q[3];
rz(-0.23770556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.33092734) q[2];
sx q[2];
rz(-1.8570447) q[2];
sx q[2];
rz(-2.1539099) q[2];
rz(0.87743131) q[3];
sx q[3];
rz(-2.8743447) q[3];
sx q[3];
rz(-2.4115244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9577318) q[0];
sx q[0];
rz(-1.6022302) q[0];
sx q[0];
rz(-3.0027332) q[0];
rz(1.214437) q[1];
sx q[1];
rz(-1.6338394) q[1];
sx q[1];
rz(-2.3249998) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3432085) q[0];
sx q[0];
rz(-1.3460165) q[0];
sx q[0];
rz(-0.019545743) q[0];
rz(2.3733489) q[2];
sx q[2];
rz(-1.381594) q[2];
sx q[2];
rz(0.18174325) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.066799435) q[1];
sx q[1];
rz(-1.5333251) q[1];
sx q[1];
rz(0.85185677) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2927033) q[3];
sx q[3];
rz(-1.0378285) q[3];
sx q[3];
rz(0.70895586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3198118) q[2];
sx q[2];
rz(-0.37415162) q[2];
sx q[2];
rz(0.4846586) q[2];
rz(-0.36568493) q[3];
sx q[3];
rz(-1.3644783) q[3];
sx q[3];
rz(-1.9827838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1300238) q[0];
sx q[0];
rz(-2.8304709) q[0];
sx q[0];
rz(3.044361) q[0];
rz(0.10487996) q[1];
sx q[1];
rz(-2.361894) q[1];
sx q[1];
rz(1.3146776) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.275248) q[0];
sx q[0];
rz(-1.5731678) q[0];
sx q[0];
rz(1.5644685) q[0];
rz(-pi) q[1];
rz(0.22366053) q[2];
sx q[2];
rz(-1.3350492) q[2];
sx q[2];
rz(-3.0652114) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7200249) q[1];
sx q[1];
rz(-2.0889258) q[1];
sx q[1];
rz(1.6314477) q[1];
rz(0.087160067) q[3];
sx q[3];
rz(-1.0356552) q[3];
sx q[3];
rz(0.59628192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.17681992) q[2];
sx q[2];
rz(-1.7440045) q[2];
sx q[2];
rz(2.9883265) q[2];
rz(-1.2166474) q[3];
sx q[3];
rz(-2.9235268) q[3];
sx q[3];
rz(-2.5254068) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1185054) q[0];
sx q[0];
rz(-0.0054792976) q[0];
sx q[0];
rz(-1.6368921) q[0];
rz(-2.3432689) q[1];
sx q[1];
rz(-1.6183805) q[1];
sx q[1];
rz(-0.33531478) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17851795) q[0];
sx q[0];
rz(-1.3097242) q[0];
sx q[0];
rz(0.065351323) q[0];
rz(1.5952571) q[2];
sx q[2];
rz(-2.3210397) q[2];
sx q[2];
rz(0.94708196) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0287231) q[1];
sx q[1];
rz(-1.2838893) q[1];
sx q[1];
rz(-0.35526459) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9934523) q[3];
sx q[3];
rz(-1.4386739) q[3];
sx q[3];
rz(-2.7293049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.12751427) q[2];
sx q[2];
rz(-1.1684343) q[2];
sx q[2];
rz(0.43759313) q[2];
rz(-1.0420927) q[3];
sx q[3];
rz(-1.7362678) q[3];
sx q[3];
rz(-0.54245814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1522778) q[0];
sx q[0];
rz(-2.2570026) q[0];
sx q[0];
rz(0.44961318) q[0];
rz(2.1647029) q[1];
sx q[1];
rz(-1.6376817) q[1];
sx q[1];
rz(-0.50122112) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4888484) q[0];
sx q[0];
rz(-2.7978659) q[0];
sx q[0];
rz(0.97596844) q[0];
x q[1];
rz(-1.1468723) q[2];
sx q[2];
rz(-1.7676643) q[2];
sx q[2];
rz(-2.9483861) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.93443643) q[1];
sx q[1];
rz(-2.0324576) q[1];
sx q[1];
rz(-2.8769879) q[1];
rz(-0.34256713) q[3];
sx q[3];
rz(-2.0122572) q[3];
sx q[3];
rz(1.6395237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0235128) q[2];
sx q[2];
rz(-2.2492275) q[2];
sx q[2];
rz(0.21772131) q[2];
rz(-1.9067541) q[3];
sx q[3];
rz(-2.4775938) q[3];
sx q[3];
rz(-0.92065221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49391838) q[0];
sx q[0];
rz(-1.7128581) q[0];
sx q[0];
rz(0.38059522) q[0];
rz(-0.71161288) q[1];
sx q[1];
rz(-1.2672392) q[1];
sx q[1];
rz(-1.738501) q[1];
rz(-3.0416476) q[2];
sx q[2];
rz(-2.6700085) q[2];
sx q[2];
rz(2.811583) q[2];
rz(-0.40583457) q[3];
sx q[3];
rz(-1.7885699) q[3];
sx q[3];
rz(-1.7287265) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
