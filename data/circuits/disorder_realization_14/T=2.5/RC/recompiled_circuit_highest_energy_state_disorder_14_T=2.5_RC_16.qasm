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
rz(2.4013588) q[0];
sx q[0];
rz(-1.6594247) q[0];
sx q[0];
rz(-2.8066714) q[0];
rz(-2.6236293) q[1];
sx q[1];
rz(-2.1393175) q[1];
sx q[1];
rz(2.5340773) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1231287) q[0];
sx q[0];
rz(-2.0573318) q[0];
sx q[0];
rz(2.8790265) q[0];
rz(0.11694853) q[2];
sx q[2];
rz(-2.0772572) q[2];
sx q[2];
rz(0.036265515) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7360037) q[1];
sx q[1];
rz(-2.0479255) q[1];
sx q[1];
rz(-1.7504343) q[1];
rz(-2.7685952) q[3];
sx q[3];
rz(-1.1524876) q[3];
sx q[3];
rz(1.4137063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.5241549) q[2];
sx q[2];
rz(-1.2158771) q[2];
sx q[2];
rz(2.4784135) q[2];
rz(-3.0607306) q[3];
sx q[3];
rz(-0.20564779) q[3];
sx q[3];
rz(1.9614356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8993768) q[0];
sx q[0];
rz(-1.3726534) q[0];
sx q[0];
rz(-0.85897613) q[0];
rz(1.8513177) q[1];
sx q[1];
rz(-1.4651508) q[1];
sx q[1];
rz(1.4069517) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55757421) q[0];
sx q[0];
rz(-1.0913335) q[0];
sx q[0];
rz(-0.25734253) q[0];
rz(-0.036769899) q[2];
sx q[2];
rz(-1.2388133) q[2];
sx q[2];
rz(0.22184243) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3000112) q[1];
sx q[1];
rz(-1.661013) q[1];
sx q[1];
rz(2.0848227) q[1];
x q[2];
rz(-2.5646016) q[3];
sx q[3];
rz(-1.9245879) q[3];
sx q[3];
rz(2.4959223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0824288) q[2];
sx q[2];
rz(-2.6291206) q[2];
sx q[2];
rz(1.814369) q[2];
rz(1.0446769) q[3];
sx q[3];
rz(-2.4278214) q[3];
sx q[3];
rz(-0.41675848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5752207) q[0];
sx q[0];
rz(-2.2307668) q[0];
sx q[0];
rz(0.4253934) q[0];
rz(-1.7644024) q[1];
sx q[1];
rz(-1.6433989) q[1];
sx q[1];
rz(-1.4345217) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4664513) q[0];
sx q[0];
rz(-1.688862) q[0];
sx q[0];
rz(2.4982014) q[0];
rz(-0.038952053) q[2];
sx q[2];
rz(-0.70868451) q[2];
sx q[2];
rz(1.5295636) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.94043865) q[1];
sx q[1];
rz(-0.80440264) q[1];
sx q[1];
rz(0.589358) q[1];
rz(-1.1528963) q[3];
sx q[3];
rz(-0.66158453) q[3];
sx q[3];
rz(-1.0583056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7003358) q[2];
sx q[2];
rz(-1.8852899) q[2];
sx q[2];
rz(-1.2909935) q[2];
rz(-1.357632) q[3];
sx q[3];
rz(-1.5057526) q[3];
sx q[3];
rz(1.4972081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5212379) q[0];
sx q[0];
rz(-1.0076032) q[0];
sx q[0];
rz(-0.77019101) q[0];
rz(2.1417446) q[1];
sx q[1];
rz(-2.5412173) q[1];
sx q[1];
rz(1.8005449) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5504549) q[0];
sx q[0];
rz(-1.7831752) q[0];
sx q[0];
rz(2.5881833) q[0];
rz(-pi) q[1];
rz(-0.88444986) q[2];
sx q[2];
rz(-1.2053066) q[2];
sx q[2];
rz(2.3523503) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8902258) q[1];
sx q[1];
rz(-2.736675) q[1];
sx q[1];
rz(-0.7271073) q[1];
rz(-pi) q[2];
x q[2];
rz(0.099319066) q[3];
sx q[3];
rz(-1.2867974) q[3];
sx q[3];
rz(-1.2885531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8009214) q[2];
sx q[2];
rz(-1.069671) q[2];
sx q[2];
rz(-0.7684024) q[2];
rz(0.10041222) q[3];
sx q[3];
rz(-1.5407591) q[3];
sx q[3];
rz(-1.8628619) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1860344) q[0];
sx q[0];
rz(-2.8362507) q[0];
sx q[0];
rz(-2.9146063) q[0];
rz(1.3817894) q[1];
sx q[1];
rz(-0.58265668) q[1];
sx q[1];
rz(1.3963799) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5433301) q[0];
sx q[0];
rz(-2.6285183) q[0];
sx q[0];
rz(0.76900478) q[0];
x q[1];
rz(-1.7368083) q[2];
sx q[2];
rz(-1.976227) q[2];
sx q[2];
rz(0.16399461) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.13670838) q[1];
sx q[1];
rz(-1.6423823) q[1];
sx q[1];
rz(-0.63386713) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.030524039) q[3];
sx q[3];
rz(-0.82477335) q[3];
sx q[3];
rz(-1.8061226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9153626) q[2];
sx q[2];
rz(-1.1504983) q[2];
sx q[2];
rz(3.0653595) q[2];
rz(-1.7539615) q[3];
sx q[3];
rz(-2.1662655) q[3];
sx q[3];
rz(1.4001728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6607894) q[0];
sx q[0];
rz(-1.5810409) q[0];
sx q[0];
rz(-1.8060818) q[0];
rz(-0.90323365) q[1];
sx q[1];
rz(-1.3390373) q[1];
sx q[1];
rz(0.18641557) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.37019) q[0];
sx q[0];
rz(-2.0982842) q[0];
sx q[0];
rz(0.50748388) q[0];
rz(-pi) q[1];
rz(-1.769563) q[2];
sx q[2];
rz(-1.0093401) q[2];
sx q[2];
rz(1.1793062) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.820436) q[1];
sx q[1];
rz(-1.3296491) q[1];
sx q[1];
rz(2.2187869) q[1];
rz(-pi) q[2];
rz(2.647688) q[3];
sx q[3];
rz(-2.0155689) q[3];
sx q[3];
rz(-0.95361036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.67941252) q[2];
sx q[2];
rz(-1.9483515) q[2];
sx q[2];
rz(-1.3235486) q[2];
rz(-1.1421674) q[3];
sx q[3];
rz(-2.3565632) q[3];
sx q[3];
rz(-0.89046684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(2.6776176) q[0];
sx q[0];
rz(-2.9587726) q[0];
sx q[0];
rz(2.5073945) q[0];
rz(1.0448666) q[1];
sx q[1];
rz(-0.8909145) q[1];
sx q[1];
rz(-2.1814836) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8329937) q[0];
sx q[0];
rz(-1.9694298) q[0];
sx q[0];
rz(-3.1118168) q[0];
rz(0.025770806) q[2];
sx q[2];
rz(-0.30116815) q[2];
sx q[2];
rz(-2.0370551) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3030678) q[1];
sx q[1];
rz(-2.6160598) q[1];
sx q[1];
rz(-1.5938894) q[1];
x q[2];
rz(2.2957357) q[3];
sx q[3];
rz(-2.1921033) q[3];
sx q[3];
rz(-1.3267856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1977957) q[2];
sx q[2];
rz(-2.4891977) q[2];
sx q[2];
rz(1.5956399) q[2];
rz(1.4173077) q[3];
sx q[3];
rz(-0.82481074) q[3];
sx q[3];
rz(-1.5203169) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5328131) q[0];
sx q[0];
rz(-2.9488035) q[0];
sx q[0];
rz(0.97484318) q[0];
rz(-0.10803647) q[1];
sx q[1];
rz(-1.8861176) q[1];
sx q[1];
rz(-1.9727762) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0625789) q[0];
sx q[0];
rz(-2.4860408) q[0];
sx q[0];
rz(2.2961839) q[0];
x q[1];
rz(1.5708099) q[2];
sx q[2];
rz(-2.332649) q[2];
sx q[2];
rz(2.1894313) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4296234) q[1];
sx q[1];
rz(-1.6584466) q[1];
sx q[1];
rz(-1.3369249) q[1];
rz(-pi) q[2];
rz(-1.5904558) q[3];
sx q[3];
rz(-1.0791313) q[3];
sx q[3];
rz(1.2439072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.03269) q[2];
sx q[2];
rz(-0.87778512) q[2];
sx q[2];
rz(-0.6558134) q[2];
rz(-0.10410318) q[3];
sx q[3];
rz(-1.3452353) q[3];
sx q[3];
rz(2.9534705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.8769237) q[0];
sx q[0];
rz(-1.3035362) q[0];
sx q[0];
rz(-1.4917829) q[0];
rz(1.0393633) q[1];
sx q[1];
rz(-2.3014258) q[1];
sx q[1];
rz(0.46844354) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7130947) q[0];
sx q[0];
rz(-1.5665496) q[0];
sx q[0];
rz(0.18761401) q[0];
rz(-pi) q[1];
rz(-1.2861112) q[2];
sx q[2];
rz(-1.9980901) q[2];
sx q[2];
rz(2.9197249) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8711618) q[1];
sx q[1];
rz(-0.79177815) q[1];
sx q[1];
rz(1.3429848) q[1];
rz(1.8519206) q[3];
sx q[3];
rz(-1.8418962) q[3];
sx q[3];
rz(0.93096126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.046772) q[2];
sx q[2];
rz(-1.5609317) q[2];
sx q[2];
rz(2.2136733) q[2];
rz(-1.3377442) q[3];
sx q[3];
rz(-2.6313621) q[3];
sx q[3];
rz(1.4891589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054166404) q[0];
sx q[0];
rz(-1.0324284) q[0];
sx q[0];
rz(-1.077865) q[0];
rz(-0.36733356) q[1];
sx q[1];
rz(-1.3307738) q[1];
sx q[1];
rz(-1.2841388) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36791641) q[0];
sx q[0];
rz(-1.8225637) q[0];
sx q[0];
rz(1.0836224) q[0];
x q[1];
rz(3.025829) q[2];
sx q[2];
rz(-2.2211005) q[2];
sx q[2];
rz(2.8011326) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.974677) q[1];
sx q[1];
rz(-1.484718) q[1];
sx q[1];
rz(-2.9831216) q[1];
x q[2];
rz(0.23350164) q[3];
sx q[3];
rz(-0.74849183) q[3];
sx q[3];
rz(-0.12602636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.0014570634) q[2];
sx q[2];
rz(-2.6900901) q[2];
sx q[2];
rz(-0.4293116) q[2];
rz(-2.9863827) q[3];
sx q[3];
rz(-2.8838172) q[3];
sx q[3];
rz(-1.3252873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24406381) q[0];
sx q[0];
rz(-1.5915992) q[0];
sx q[0];
rz(-1.5912548) q[0];
rz(-2.2776729) q[1];
sx q[1];
rz(-2.7680631) q[1];
sx q[1];
rz(1.6815129) q[1];
rz(2.7511394) q[2];
sx q[2];
rz(-0.57885546) q[2];
sx q[2];
rz(-0.59138966) q[2];
rz(-0.33959099) q[3];
sx q[3];
rz(-0.97151269) q[3];
sx q[3];
rz(1.1480939) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
