OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.03117938) q[0];
sx q[0];
rz(-0.92209446) q[0];
sx q[0];
rz(2.3715012) q[0];
rz(0.71647477) q[1];
sx q[1];
rz(-2.1991576) q[1];
sx q[1];
rz(0.64185774) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7866369) q[0];
sx q[0];
rz(-1.6842951) q[0];
sx q[0];
rz(0.20239511) q[0];
x q[1];
rz(1.5719169) q[2];
sx q[2];
rz(-1.5693451) q[2];
sx q[2];
rz(3.0641132) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.042965021) q[1];
sx q[1];
rz(-2.7471099) q[1];
sx q[1];
rz(-2.7261655) q[1];
rz(2.558488) q[3];
sx q[3];
rz(-0.92188406) q[3];
sx q[3];
rz(-0.16498868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1610819) q[2];
sx q[2];
rz(-1.0970205) q[2];
sx q[2];
rz(-0.80992997) q[2];
rz(0.03446456) q[3];
sx q[3];
rz(-2.4813215) q[3];
sx q[3];
rz(-3.0380429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63475364) q[0];
sx q[0];
rz(-0.50951183) q[0];
sx q[0];
rz(0.65938812) q[0];
rz(-1.7022645) q[1];
sx q[1];
rz(-1.5123475) q[1];
sx q[1];
rz(2.6606681) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60431193) q[0];
sx q[0];
rz(-1.7780281) q[0];
sx q[0];
rz(-2.0181433) q[0];
x q[1];
rz(2.2873053) q[2];
sx q[2];
rz(-2.169974) q[2];
sx q[2];
rz(-2.6024659) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9316599) q[1];
sx q[1];
rz(-1.6012962) q[1];
sx q[1];
rz(-1.362544) q[1];
rz(-pi) q[2];
rz(-1.102785) q[3];
sx q[3];
rz(-1.0429405) q[3];
sx q[3];
rz(-0.46609391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3769569) q[2];
sx q[2];
rz(-0.83676338) q[2];
sx q[2];
rz(-2.5118206) q[2];
rz(1.9624286) q[3];
sx q[3];
rz(-0.7032913) q[3];
sx q[3];
rz(-1.111697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69029194) q[0];
sx q[0];
rz(-3.1307104) q[0];
sx q[0];
rz(-0.96845281) q[0];
rz(0.14006607) q[1];
sx q[1];
rz(-1.3532956) q[1];
sx q[1];
rz(0.5828988) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.934865) q[0];
sx q[0];
rz(-1.7375653) q[0];
sx q[0];
rz(-1.5251446) q[0];
x q[1];
rz(0.83816041) q[2];
sx q[2];
rz(-1.7485022) q[2];
sx q[2];
rz(-0.017824307) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0839094) q[1];
sx q[1];
rz(-2.1321802) q[1];
sx q[1];
rz(-2.3777005) q[1];
x q[2];
rz(-0.85286136) q[3];
sx q[3];
rz(-0.4133458) q[3];
sx q[3];
rz(1.0173544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.53050238) q[2];
sx q[2];
rz(-0.077433057) q[2];
sx q[2];
rz(0.21406847) q[2];
rz(0.30238447) q[3];
sx q[3];
rz(-2.3979135) q[3];
sx q[3];
rz(-2.7157057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9532303) q[0];
sx q[0];
rz(-2.132405) q[0];
sx q[0];
rz(2.0971712) q[0];
rz(-2.2638679) q[1];
sx q[1];
rz(-1.5889771) q[1];
sx q[1];
rz(0.16876076) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.343086) q[0];
sx q[0];
rz(-2.0650535) q[0];
sx q[0];
rz(-2.8841208) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8977106) q[2];
sx q[2];
rz(-1.765328) q[2];
sx q[2];
rz(-3.1049797) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5280594) q[1];
sx q[1];
rz(-2.4984887) q[1];
sx q[1];
rz(-0.36524857) q[1];
rz(-pi) q[2];
rz(3.0996347) q[3];
sx q[3];
rz(-2.2760512) q[3];
sx q[3];
rz(2.369997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7793444) q[2];
sx q[2];
rz(-1.2835953) q[2];
sx q[2];
rz(-2.0812422) q[2];
rz(0.29252163) q[3];
sx q[3];
rz(-0.7258324) q[3];
sx q[3];
rz(-0.084107548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5523858) q[0];
sx q[0];
rz(-0.64738208) q[0];
sx q[0];
rz(2.2733083) q[0];
rz(2.5153416) q[1];
sx q[1];
rz(-1.3887082) q[1];
sx q[1];
rz(-1.7832322) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6407961) q[0];
sx q[0];
rz(-1.2828886) q[0];
sx q[0];
rz(-1.663371) q[0];
rz(2.7600438) q[2];
sx q[2];
rz(-1.7722436) q[2];
sx q[2];
rz(0.82009456) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7029801) q[1];
sx q[1];
rz(-1.3452521) q[1];
sx q[1];
rz(1.799052) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63774469) q[3];
sx q[3];
rz(-1.2079835) q[3];
sx q[3];
rz(2.2408298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2061578) q[2];
sx q[2];
rz(-0.81659603) q[2];
sx q[2];
rz(2.6182776) q[2];
rz(-2.7715136) q[3];
sx q[3];
rz(-0.78032929) q[3];
sx q[3];
rz(-1.3453329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0321781) q[0];
sx q[0];
rz(-2.9258756) q[0];
sx q[0];
rz(-0.35414645) q[0];
rz(2.1996563) q[1];
sx q[1];
rz(-1.4553921) q[1];
sx q[1];
rz(-1.8249493) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9958651) q[0];
sx q[0];
rz(-0.52477974) q[0];
sx q[0];
rz(-0.20916931) q[0];
rz(-pi) q[1];
rz(1.5572335) q[2];
sx q[2];
rz(-2.1199653) q[2];
sx q[2];
rz(-1.4209335) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5806953) q[1];
sx q[1];
rz(-1.744701) q[1];
sx q[1];
rz(-1.5087125) q[1];
rz(-pi) q[2];
rz(1.2888258) q[3];
sx q[3];
rz(-2.1957955) q[3];
sx q[3];
rz(0.19240141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3256623) q[2];
sx q[2];
rz(-0.38620913) q[2];
sx q[2];
rz(2.8032934) q[2];
rz(-0.47652388) q[3];
sx q[3];
rz(-2.3881113) q[3];
sx q[3];
rz(2.8009955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7020096) q[0];
sx q[0];
rz(-1.5832573) q[0];
sx q[0];
rz(2.6038792) q[0];
rz(1.4078377) q[1];
sx q[1];
rz(-0.45999637) q[1];
sx q[1];
rz(-0.62526155) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38442507) q[0];
sx q[0];
rz(-1.0119857) q[0];
sx q[0];
rz(1.998087) q[0];
rz(-pi) q[1];
rz(1.4971259) q[2];
sx q[2];
rz(-2.1084614) q[2];
sx q[2];
rz(-1.0755838) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3846674) q[1];
sx q[1];
rz(-2.1672492) q[1];
sx q[1];
rz(-1.8310962) q[1];
x q[2];
rz(-1.2042562) q[3];
sx q[3];
rz(-0.82895422) q[3];
sx q[3];
rz(-0.77223611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9739146) q[2];
sx q[2];
rz(-2.6748071) q[2];
sx q[2];
rz(-1.7612877) q[2];
rz(2.7050833) q[3];
sx q[3];
rz(-2.1197539) q[3];
sx q[3];
rz(2.4280587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1768782) q[0];
sx q[0];
rz(-2.5202993) q[0];
sx q[0];
rz(-2.6253413) q[0];
rz(-2.3618354) q[1];
sx q[1];
rz(-2.1596491) q[1];
sx q[1];
rz(-2.0947184) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6927495) q[0];
sx q[0];
rz(-1.4268645) q[0];
sx q[0];
rz(2.8400665) q[0];
x q[1];
rz(2.4372836) q[2];
sx q[2];
rz(-0.71137911) q[2];
sx q[2];
rz(-0.34397438) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2566764) q[1];
sx q[1];
rz(-3.0312928) q[1];
sx q[1];
rz(-1.0156471) q[1];
x q[2];
rz(-0.67892142) q[3];
sx q[3];
rz(-2.553294) q[3];
sx q[3];
rz(-3.1117518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9376935) q[2];
sx q[2];
rz(-2.9475309) q[2];
sx q[2];
rz(0.95721179) q[2];
rz(-2.8474478) q[3];
sx q[3];
rz(-1.1295986) q[3];
sx q[3];
rz(0.2778151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1421563) q[0];
sx q[0];
rz(-0.23040982) q[0];
sx q[0];
rz(0.18705046) q[0];
rz(2.710178) q[1];
sx q[1];
rz(-2.7038733) q[1];
sx q[1];
rz(1.6492708) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87875116) q[0];
sx q[0];
rz(-3.0410112) q[0];
sx q[0];
rz(-0.79128964) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9052986) q[2];
sx q[2];
rz(-1.750947) q[2];
sx q[2];
rz(1.0049154) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.88384151) q[1];
sx q[1];
rz(-1.7353936) q[1];
sx q[1];
rz(1.1562892) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.047916) q[3];
sx q[3];
rz(-1.4149349) q[3];
sx q[3];
rz(0.81934281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.46818587) q[2];
sx q[2];
rz(-2.2298721) q[2];
sx q[2];
rz(2.1569596) q[2];
rz(0.51472384) q[3];
sx q[3];
rz(-0.52557164) q[3];
sx q[3];
rz(1.2334067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.89727) q[0];
sx q[0];
rz(-1.4800625) q[0];
sx q[0];
rz(-2.3829714) q[0];
rz(1.1933391) q[1];
sx q[1];
rz(-1.9494282) q[1];
sx q[1];
rz(1.4512482) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0139931) q[0];
sx q[0];
rz(-1.3422478) q[0];
sx q[0];
rz(2.9319113) q[0];
rz(-pi) q[1];
rz(-1.7733943) q[2];
sx q[2];
rz(-2.9349064) q[2];
sx q[2];
rz(0.24487534) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9368694) q[1];
sx q[1];
rz(-1.544049) q[1];
sx q[1];
rz(2.1863947) q[1];
rz(-pi) q[2];
rz(1.9791466) q[3];
sx q[3];
rz(-1.8682042) q[3];
sx q[3];
rz(1.4434467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.54185581) q[2];
sx q[2];
rz(-0.92076045) q[2];
sx q[2];
rz(2.3403781) q[2];
rz(-2.2090705) q[3];
sx q[3];
rz(-1.9263809) q[3];
sx q[3];
rz(0.41509375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5634609) q[0];
sx q[0];
rz(-1.3787855) q[0];
sx q[0];
rz(-0.61510573) q[0];
rz(2.8201132) q[1];
sx q[1];
rz(-0.97578661) q[1];
sx q[1];
rz(-1.6937561) q[1];
rz(-1.4463439) q[2];
sx q[2];
rz(-1.7348518) q[2];
sx q[2];
rz(3.0748925) q[2];
rz(-1.2899756) q[3];
sx q[3];
rz(-1.7924037) q[3];
sx q[3];
rz(-2.0879346) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
