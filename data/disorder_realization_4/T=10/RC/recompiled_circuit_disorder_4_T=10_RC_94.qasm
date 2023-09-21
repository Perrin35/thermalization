OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.714158) q[0];
sx q[0];
rz(-2.7058869) q[0];
sx q[0];
rz(-0.92619196) q[0];
rz(1.9594833) q[1];
sx q[1];
rz(-0.73298454) q[1];
sx q[1];
rz(-2.7690673) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0988136) q[0];
sx q[0];
rz(-0.61646898) q[0];
sx q[0];
rz(-0.35538506) q[0];
rz(-pi) q[1];
rz(-1.4650605) q[2];
sx q[2];
rz(-0.9423965) q[2];
sx q[2];
rz(-2.489593) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8718611) q[1];
sx q[1];
rz(-1.6850123) q[1];
sx q[1];
rz(1.3807339) q[1];
rz(-pi) q[2];
rz(-0.11301179) q[3];
sx q[3];
rz(-1.3556619) q[3];
sx q[3];
rz(-2.5290031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.42840502) q[2];
sx q[2];
rz(-1.5593854) q[2];
sx q[2];
rz(-0.92450809) q[2];
rz(1.472578) q[3];
sx q[3];
rz(-1.893483) q[3];
sx q[3];
rz(1.4991466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18838841) q[0];
sx q[0];
rz(-0.062347118) q[0];
sx q[0];
rz(-1.5361319) q[0];
rz(-2.9470782) q[1];
sx q[1];
rz(-1.3214) q[1];
sx q[1];
rz(0.054873437) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40614906) q[0];
sx q[0];
rz(-0.14649728) q[0];
sx q[0];
rz(-0.88645257) q[0];
x q[1];
rz(0.61972159) q[2];
sx q[2];
rz(-1.1148858) q[2];
sx q[2];
rz(1.0085307) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84150746) q[1];
sx q[1];
rz(-1.9133139) q[1];
sx q[1];
rz(0.0042580982) q[1];
rz(-0.26344928) q[3];
sx q[3];
rz(-2.4512495) q[3];
sx q[3];
rz(-0.024397959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.43869552) q[2];
sx q[2];
rz(-2.5800811) q[2];
sx q[2];
rz(-1.4820209) q[2];
rz(-2.1510018) q[3];
sx q[3];
rz(-1.7329268) q[3];
sx q[3];
rz(-1.7747442) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7822587) q[0];
sx q[0];
rz(-0.5193091) q[0];
sx q[0];
rz(-0.3749795) q[0];
rz(0.24770501) q[1];
sx q[1];
rz(-1.9539555) q[1];
sx q[1];
rz(-1.3365655) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35555392) q[0];
sx q[0];
rz(-0.22909129) q[0];
sx q[0];
rz(-1.3130207) q[0];
rz(-pi) q[1];
rz(-0.58972085) q[2];
sx q[2];
rz(-1.1955185) q[2];
sx q[2];
rz(2.7207665) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6429236) q[1];
sx q[1];
rz(-1.8734697) q[1];
sx q[1];
rz(-2.32294) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3744208) q[3];
sx q[3];
rz(-1.2734405) q[3];
sx q[3];
rz(1.5470099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0388564) q[2];
sx q[2];
rz(-0.83013022) q[2];
sx q[2];
rz(2.1170199) q[2];
rz(-1.864795) q[3];
sx q[3];
rz(-1.5921311) q[3];
sx q[3];
rz(-1.6259441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17722002) q[0];
sx q[0];
rz(-1.6765046) q[0];
sx q[0];
rz(3.1080416) q[0];
rz(-1.9112446) q[1];
sx q[1];
rz(-0.76554811) q[1];
sx q[1];
rz(1.7376602) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0899635) q[0];
sx q[0];
rz(-1.7647867) q[0];
sx q[0];
rz(-2.7882663) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3061182) q[2];
sx q[2];
rz(-0.88047853) q[2];
sx q[2];
rz(-0.26963216) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2851706) q[1];
sx q[1];
rz(-0.31177786) q[1];
sx q[1];
rz(-1.9429893) q[1];
rz(-pi) q[2];
rz(1.674026) q[3];
sx q[3];
rz(-0.69403115) q[3];
sx q[3];
rz(-0.3895143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4924865) q[2];
sx q[2];
rz(-2.2968569) q[2];
sx q[2];
rz(0.21326324) q[2];
rz(0.02031859) q[3];
sx q[3];
rz(-1.3213986) q[3];
sx q[3];
rz(0.4030574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3605109) q[0];
sx q[0];
rz(-1.1163982) q[0];
sx q[0];
rz(2.1257341) q[0];
rz(-1.6332743) q[1];
sx q[1];
rz(-0.95817482) q[1];
sx q[1];
rz(-3.1075081) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8071472) q[0];
sx q[0];
rz(-0.42754081) q[0];
sx q[0];
rz(-2.7464675) q[0];
rz(2.5675315) q[2];
sx q[2];
rz(-1.7957557) q[2];
sx q[2];
rz(-0.056919295) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.34921131) q[1];
sx q[1];
rz(-1.5639018) q[1];
sx q[1];
rz(-1.8551808) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91655101) q[3];
sx q[3];
rz(-0.74186462) q[3];
sx q[3];
rz(1.3548917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.74026996) q[2];
sx q[2];
rz(-0.57381845) q[2];
sx q[2];
rz(1.9704698) q[2];
rz(1.3327538) q[3];
sx q[3];
rz(-1.4274024) q[3];
sx q[3];
rz(-0.12935054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68833441) q[0];
sx q[0];
rz(-1.1880705) q[0];
sx q[0];
rz(0.43584287) q[0];
rz(-2.0360937) q[1];
sx q[1];
rz(-1.2970129) q[1];
sx q[1];
rz(1.9357392) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3656475) q[0];
sx q[0];
rz(-1.8071113) q[0];
sx q[0];
rz(-0.81415557) q[0];
x q[1];
rz(-0.016024307) q[2];
sx q[2];
rz(-1.4756225) q[2];
sx q[2];
rz(1.3118088) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7534415) q[1];
sx q[1];
rz(-2.6302164) q[1];
sx q[1];
rz(-3.0877153) q[1];
rz(-pi) q[2];
rz(0.033127012) q[3];
sx q[3];
rz(-1.1713542) q[3];
sx q[3];
rz(-1.9314194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.897052) q[2];
sx q[2];
rz(-2.5138833) q[2];
sx q[2];
rz(-0.051606027) q[2];
rz(-2.7339325) q[3];
sx q[3];
rz(-1.1838653) q[3];
sx q[3];
rz(2.6962962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62149858) q[0];
sx q[0];
rz(-0.61722732) q[0];
sx q[0];
rz(-2.4682585) q[0];
rz(-0.78397059) q[1];
sx q[1];
rz(-1.7242804) q[1];
sx q[1];
rz(-0.53692445) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4250701) q[0];
sx q[0];
rz(-1.5050097) q[0];
sx q[0];
rz(1.6545047) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0836673) q[2];
sx q[2];
rz(-0.29209902) q[2];
sx q[2];
rz(-0.7823173) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9432224) q[1];
sx q[1];
rz(-0.80723019) q[1];
sx q[1];
rz(3.0109349) q[1];
rz(-pi) q[2];
rz(0.59099205) q[3];
sx q[3];
rz(-1.4431074) q[3];
sx q[3];
rz(-1.5203116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9817104) q[2];
sx q[2];
rz(-1.8843702) q[2];
sx q[2];
rz(2.9505777) q[2];
rz(-2.8602709) q[3];
sx q[3];
rz(-1.9502935) q[3];
sx q[3];
rz(0.34240001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1087588) q[0];
sx q[0];
rz(-0.37029752) q[0];
sx q[0];
rz(-0.38761815) q[0];
rz(-0.11501137) q[1];
sx q[1];
rz(-1.7128046) q[1];
sx q[1];
rz(-0.33755916) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96466366) q[0];
sx q[0];
rz(-1.7690072) q[0];
sx q[0];
rz(-0.49074178) q[0];
x q[1];
rz(0.07143306) q[2];
sx q[2];
rz(-0.93284235) q[2];
sx q[2];
rz(-0.86779867) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.464752) q[1];
sx q[1];
rz(-2.5555579) q[1];
sx q[1];
rz(1.3321484) q[1];
x q[2];
rz(-3.0739325) q[3];
sx q[3];
rz(-0.75660556) q[3];
sx q[3];
rz(-2.9697231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9541624) q[2];
sx q[2];
rz(-0.53415161) q[2];
sx q[2];
rz(1.2672651) q[2];
rz(-0.77504843) q[3];
sx q[3];
rz(-1.5379484) q[3];
sx q[3];
rz(0.78139853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18701126) q[0];
sx q[0];
rz(-2.1033852) q[0];
sx q[0];
rz(-2.9206081) q[0];
rz(0.9221319) q[1];
sx q[1];
rz(-1.2698959) q[1];
sx q[1];
rz(-1.0029213) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.636165) q[0];
sx q[0];
rz(-1.822585) q[0];
sx q[0];
rz(-1.4279143) q[0];
rz(-2.2641364) q[2];
sx q[2];
rz(-0.71091953) q[2];
sx q[2];
rz(0.029821776) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.50833118) q[1];
sx q[1];
rz(-1.08053) q[1];
sx q[1];
rz(2.0983178) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8691611) q[3];
sx q[3];
rz(-2.2386132) q[3];
sx q[3];
rz(1.2215134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6691436) q[2];
sx q[2];
rz(-2.3995235) q[2];
sx q[2];
rz(2.9830902) q[2];
rz(2.6878099) q[3];
sx q[3];
rz(-0.78275371) q[3];
sx q[3];
rz(-2.2973072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6145988) q[0];
sx q[0];
rz(-0.90899962) q[0];
sx q[0];
rz(-2.9504543) q[0];
rz(-2.846431) q[1];
sx q[1];
rz(-2.252153) q[1];
sx q[1];
rz(2.2492762) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44360456) q[0];
sx q[0];
rz(-0.83570489) q[0];
sx q[0];
rz(0.23905917) q[0];
rz(-2.2116682) q[2];
sx q[2];
rz(-1.6460437) q[2];
sx q[2];
rz(-1.4571232) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1103507) q[1];
sx q[1];
rz(-2.894245) q[1];
sx q[1];
rz(-0.27242839) q[1];
rz(-0.51384135) q[3];
sx q[3];
rz(-2.6548879) q[3];
sx q[3];
rz(1.0443618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3698547) q[2];
sx q[2];
rz(-0.83738804) q[2];
sx q[2];
rz(-0.85956335) q[2];
rz(1.2236979) q[3];
sx q[3];
rz(-0.92646354) q[3];
sx q[3];
rz(-0.74444509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0201465) q[0];
sx q[0];
rz(-2.2059724) q[0];
sx q[0];
rz(-1.7511517) q[0];
rz(-1.4795115) q[1];
sx q[1];
rz(-2.6585487) q[1];
sx q[1];
rz(-1.9225635) q[1];
rz(-1.4798726) q[2];
sx q[2];
rz(-0.29279136) q[2];
sx q[2];
rz(1.472483) q[2];
rz(-1.0449833) q[3];
sx q[3];
rz(-3.0621739) q[3];
sx q[3];
rz(-0.0041181507) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
