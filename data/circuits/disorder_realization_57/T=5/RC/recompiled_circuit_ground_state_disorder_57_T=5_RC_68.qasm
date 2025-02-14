OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.39855555) q[0];
sx q[0];
rz(-0.86657137) q[0];
sx q[0];
rz(-2.804629) q[0];
rz(0.26951867) q[1];
sx q[1];
rz(4.0842847) q[1];
sx q[1];
rz(10.684291) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.058858697) q[0];
sx q[0];
rz(-1.0089416) q[0];
sx q[0];
rz(-2.492172) q[0];
x q[1];
rz(1.3607698) q[2];
sx q[2];
rz(-2.2968596) q[2];
sx q[2];
rz(-0.10440102) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0639021) q[1];
sx q[1];
rz(-0.55455583) q[1];
sx q[1];
rz(0.18973236) q[1];
x q[2];
rz(-0.97136949) q[3];
sx q[3];
rz(-1.9183591) q[3];
sx q[3];
rz(-2.3088648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.17244615) q[2];
sx q[2];
rz(-1.5643876) q[2];
sx q[2];
rz(-1.2338314) q[2];
rz(-1.5305758) q[3];
sx q[3];
rz(-1.4065892) q[3];
sx q[3];
rz(0.56510258) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81561404) q[0];
sx q[0];
rz(-2.1765206) q[0];
sx q[0];
rz(2.8976029) q[0];
rz(-1.9832393) q[1];
sx q[1];
rz(-2.127425) q[1];
sx q[1];
rz(0.44833952) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8566003) q[0];
sx q[0];
rz(-1.0831079) q[0];
sx q[0];
rz(1.7708805) q[0];
rz(-pi) q[1];
rz(0.59780663) q[2];
sx q[2];
rz(-2.2181244) q[2];
sx q[2];
rz(1.7773599) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.068591) q[1];
sx q[1];
rz(-1.0319971) q[1];
sx q[1];
rz(-2.6122142) q[1];
rz(-pi) q[2];
rz(0.86037029) q[3];
sx q[3];
rz(-0.67796889) q[3];
sx q[3];
rz(1.9063661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9490732) q[2];
sx q[2];
rz(-3.101888) q[2];
sx q[2];
rz(-2.330759) q[2];
rz(3.132498) q[3];
sx q[3];
rz(-1.5261212) q[3];
sx q[3];
rz(0.31753376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72967616) q[0];
sx q[0];
rz(-1.7971973) q[0];
sx q[0];
rz(0.94773951) q[0];
rz(-2.2432227) q[1];
sx q[1];
rz(-2.4285474) q[1];
sx q[1];
rz(2.9352442) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8074051) q[0];
sx q[0];
rz(-1.5631521) q[0];
sx q[0];
rz(-0.81595465) q[0];
rz(-pi) q[1];
rz(-0.79282615) q[2];
sx q[2];
rz(-1.2812876) q[2];
sx q[2];
rz(-0.0057366554) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.84103407) q[1];
sx q[1];
rz(-1.3649458) q[1];
sx q[1];
rz(1.1080075) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.523866) q[3];
sx q[3];
rz(-1.8968256) q[3];
sx q[3];
rz(0.24229953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.72948939) q[2];
sx q[2];
rz(-2.1495843) q[2];
sx q[2];
rz(0.98423973) q[2];
rz(2.6442243) q[3];
sx q[3];
rz(-1.2593185) q[3];
sx q[3];
rz(2.3965912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3618149) q[0];
sx q[0];
rz(-3.0497157) q[0];
sx q[0];
rz(-0.20633695) q[0];
rz(1.4641209) q[1];
sx q[1];
rz(-2.3787777) q[1];
sx q[1];
rz(0.62613553) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7595644) q[0];
sx q[0];
rz(-2.3256178) q[0];
sx q[0];
rz(2.0618646) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5690313) q[2];
sx q[2];
rz(-1.4680099) q[2];
sx q[2];
rz(2.4461164) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4724222) q[1];
sx q[1];
rz(-1.0423933) q[1];
sx q[1];
rz(1.0086568) q[1];
x q[2];
rz(1.3317165) q[3];
sx q[3];
rz(-2.7446973) q[3];
sx q[3];
rz(-1.797685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3760486) q[2];
sx q[2];
rz(-0.77771336) q[2];
sx q[2];
rz(-0.98975873) q[2];
rz(1.0342213) q[3];
sx q[3];
rz(-1.1690305) q[3];
sx q[3];
rz(2.6160252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.702221) q[0];
sx q[0];
rz(-1.419743) q[0];
sx q[0];
rz(-1.1628344) q[0];
rz(-2.974466) q[1];
sx q[1];
rz(-1.5252599) q[1];
sx q[1];
rz(-2.4662245) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5655095) q[0];
sx q[0];
rz(-1.7890437) q[0];
sx q[0];
rz(-2.4714787) q[0];
rz(-2.8020892) q[2];
sx q[2];
rz(-1.9304747) q[2];
sx q[2];
rz(2.8946257) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.31084222) q[1];
sx q[1];
rz(-1.1427587) q[1];
sx q[1];
rz(-0.38106783) q[1];
x q[2];
rz(2.1618103) q[3];
sx q[3];
rz(-1.3434935) q[3];
sx q[3];
rz(2.1066372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2165788) q[2];
sx q[2];
rz(-1.4224956) q[2];
sx q[2];
rz(0.50755802) q[2];
rz(-0.01928586) q[3];
sx q[3];
rz(-0.79500335) q[3];
sx q[3];
rz(1.219334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60292202) q[0];
sx q[0];
rz(-2.5141073) q[0];
sx q[0];
rz(2.4399309) q[0];
rz(2.498846) q[1];
sx q[1];
rz(-1.1593436) q[1];
sx q[1];
rz(-1.8205551) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3747201) q[0];
sx q[0];
rz(-1.5581616) q[0];
sx q[0];
rz(-2.1983474) q[0];
rz(0.8708839) q[2];
sx q[2];
rz(-1.8485763) q[2];
sx q[2];
rz(-0.50123614) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.031547) q[1];
sx q[1];
rz(-2.4737691) q[1];
sx q[1];
rz(0.948443) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8529529) q[3];
sx q[3];
rz(-2.309679) q[3];
sx q[3];
rz(0.4543685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8464437) q[2];
sx q[2];
rz(-2.6602793) q[2];
sx q[2];
rz(-0.63977891) q[2];
rz(-2.4844737) q[3];
sx q[3];
rz(-1.492447) q[3];
sx q[3];
rz(2.8583728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.101864) q[0];
sx q[0];
rz(-1.855408) q[0];
sx q[0];
rz(-2.3611327) q[0];
rz(1.9038433) q[1];
sx q[1];
rz(-1.8433808) q[1];
sx q[1];
rz(0.16403988) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6934768) q[0];
sx q[0];
rz(-2.0882029) q[0];
sx q[0];
rz(0.2736926) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1272993) q[2];
sx q[2];
rz(-1.9464543) q[2];
sx q[2];
rz(-1.2835447) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.32719993) q[1];
sx q[1];
rz(-2.3770824) q[1];
sx q[1];
rz(2.3491377) q[1];
rz(2.6884671) q[3];
sx q[3];
rz(-2.2966566) q[3];
sx q[3];
rz(0.22900861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8728472) q[2];
sx q[2];
rz(-1.9445323) q[2];
sx q[2];
rz(-0.8806814) q[2];
rz(-1.667048) q[3];
sx q[3];
rz(-1.0677974) q[3];
sx q[3];
rz(-1.7269945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47377652) q[0];
sx q[0];
rz(-2.9528463) q[0];
sx q[0];
rz(2.012398) q[0];
rz(0.95343626) q[1];
sx q[1];
rz(-0.9402746) q[1];
sx q[1];
rz(-0.58852351) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1746219) q[0];
sx q[0];
rz(-1.0292639) q[0];
sx q[0];
rz(-2.7906899) q[0];
x q[1];
rz(-1.2382201) q[2];
sx q[2];
rz(-1.3540478) q[2];
sx q[2];
rz(0.34453604) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2413832) q[1];
sx q[1];
rz(-1.2523249) q[1];
sx q[1];
rz(-2.4584944) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7950012) q[3];
sx q[3];
rz(-1.5191169) q[3];
sx q[3];
rz(1.3286852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.15527655) q[2];
sx q[2];
rz(-2.1337815) q[2];
sx q[2];
rz(-2.0493832) q[2];
rz(0.80882597) q[3];
sx q[3];
rz(-1.0816962) q[3];
sx q[3];
rz(-1.4441747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5876193) q[0];
sx q[0];
rz(-1.076979) q[0];
sx q[0];
rz(2.4699566) q[0];
rz(-0.49297586) q[1];
sx q[1];
rz(-2.2608345) q[1];
sx q[1];
rz(-1.902045) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17104471) q[0];
sx q[0];
rz(-1.8080304) q[0];
sx q[0];
rz(0.0026654442) q[0];
rz(-0.19723326) q[2];
sx q[2];
rz(-1.5805564) q[2];
sx q[2];
rz(-0.41744864) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.904027) q[1];
sx q[1];
rz(-0.75566245) q[1];
sx q[1];
rz(2.4811238) q[1];
x q[2];
rz(2.3726241) q[3];
sx q[3];
rz(-1.4533224) q[3];
sx q[3];
rz(-0.75107915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.12941831) q[2];
sx q[2];
rz(-1.47379) q[2];
sx q[2];
rz(-0.54656452) q[2];
rz(2.2187345) q[3];
sx q[3];
rz(-0.62896362) q[3];
sx q[3];
rz(-1.1465237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4839812) q[0];
sx q[0];
rz(-2.681499) q[0];
sx q[0];
rz(-0.42903236) q[0];
rz(-1.8589164) q[1];
sx q[1];
rz(-1.7762643) q[1];
sx q[1];
rz(0.081258953) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6264296) q[0];
sx q[0];
rz(-1.7395955) q[0];
sx q[0];
rz(1.4909046) q[0];
rz(-pi) q[1];
rz(1.6363898) q[2];
sx q[2];
rz(-0.4609522) q[2];
sx q[2];
rz(0.84921244) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.024549896) q[1];
sx q[1];
rz(-2.3094842) q[1];
sx q[1];
rz(-2.1348272) q[1];
rz(-pi) q[2];
rz(-0.83286442) q[3];
sx q[3];
rz(-2.500415) q[3];
sx q[3];
rz(-2.8744551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.73605865) q[2];
sx q[2];
rz(-2.2106876) q[2];
sx q[2];
rz(1.1116568) q[2];
rz(-1.9723655) q[3];
sx q[3];
rz(-2.4284095) q[3];
sx q[3];
rz(-2.9986103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1375167) q[0];
sx q[0];
rz(-0.62340323) q[0];
sx q[0];
rz(2.2829983) q[0];
rz(0.89523347) q[1];
sx q[1];
rz(-1.3139071) q[1];
sx q[1];
rz(0.039176686) q[1];
rz(3.1298755) q[2];
sx q[2];
rz(-1.313971) q[2];
sx q[2];
rz(1.32709) q[2];
rz(0.17651996) q[3];
sx q[3];
rz(-0.64809496) q[3];
sx q[3];
rz(2.5213836) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
