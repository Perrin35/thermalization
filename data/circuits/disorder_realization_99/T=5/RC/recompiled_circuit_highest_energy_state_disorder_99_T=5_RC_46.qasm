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
rz(1.4003657) q[0];
sx q[0];
rz(-0.23973149) q[0];
sx q[0];
rz(0.32455197) q[0];
rz(5.3311081) q[1];
sx q[1];
rz(6.0573112) q[1];
sx q[1];
rz(4.4499302) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0265147) q[0];
sx q[0];
rz(-0.9008207) q[0];
sx q[0];
rz(-1.3967394) q[0];
x q[1];
rz(-0.24271528) q[2];
sx q[2];
rz(-2.6449892) q[2];
sx q[2];
rz(0.14641031) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.54656279) q[1];
sx q[1];
rz(-1.1076704) q[1];
sx q[1];
rz(-0.28966499) q[1];
x q[2];
rz(-2.4423787) q[3];
sx q[3];
rz(-1.4627856) q[3];
sx q[3];
rz(-2.8349769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41210458) q[2];
sx q[2];
rz(-1.8821913) q[2];
sx q[2];
rz(-0.35120249) q[2];
rz(1.8515733) q[3];
sx q[3];
rz(-1.131564) q[3];
sx q[3];
rz(1.2235519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38607645) q[0];
sx q[0];
rz(-1.2048683) q[0];
sx q[0];
rz(-2.2112041) q[0];
rz(1.3049841) q[1];
sx q[1];
rz(-0.84140673) q[1];
sx q[1];
rz(-2.5321541) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064548858) q[0];
sx q[0];
rz(-1.4676369) q[0];
sx q[0];
rz(1.9984238) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8687042) q[2];
sx q[2];
rz(-1.6570083) q[2];
sx q[2];
rz(2.6484368) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.37928615) q[1];
sx q[1];
rz(-1.2067144) q[1];
sx q[1];
rz(1.9696139) q[1];
rz(-pi) q[2];
rz(1.8279161) q[3];
sx q[3];
rz(-0.58565307) q[3];
sx q[3];
rz(-2.3975092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.096752) q[2];
sx q[2];
rz(-0.98016206) q[2];
sx q[2];
rz(3.0668104) q[2];
rz(2.6212202) q[3];
sx q[3];
rz(-2.4986391) q[3];
sx q[3];
rz(2.5274932) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5388913) q[0];
sx q[0];
rz(-0.7170054) q[0];
sx q[0];
rz(-0.30211788) q[0];
rz(0.27613861) q[1];
sx q[1];
rz(-1.3172251) q[1];
sx q[1];
rz(-1.4322697) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0574422) q[0];
sx q[0];
rz(-2.1085116) q[0];
sx q[0];
rz(2.6938963) q[0];
rz(-1.8319857) q[2];
sx q[2];
rz(-2.1495594) q[2];
sx q[2];
rz(-2.8566993) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.773743) q[1];
sx q[1];
rz(-2.5605695) q[1];
sx q[1];
rz(1.8598721) q[1];
x q[2];
rz(-3.0705419) q[3];
sx q[3];
rz(-1.1059575) q[3];
sx q[3];
rz(2.0185061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3064731) q[2];
sx q[2];
rz(-2.6758631) q[2];
sx q[2];
rz(-1.2960557) q[2];
rz(-0.45201388) q[3];
sx q[3];
rz(-1.1072423) q[3];
sx q[3];
rz(-0.38255102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6119659) q[0];
sx q[0];
rz(-1.9744248) q[0];
sx q[0];
rz(-1.8091328) q[0];
rz(-1.0491071) q[1];
sx q[1];
rz(-1.0915979) q[1];
sx q[1];
rz(1.5886935) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5063254) q[0];
sx q[0];
rz(-2.679279) q[0];
sx q[0];
rz(-0.72700951) q[0];
rz(-pi) q[1];
rz(-1.5459486) q[2];
sx q[2];
rz(-2.6406807) q[2];
sx q[2];
rz(-0.76297578) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5415948) q[1];
sx q[1];
rz(-1.3728956) q[1];
sx q[1];
rz(-2.2210414) q[1];
x q[2];
rz(2.6060206) q[3];
sx q[3];
rz(-2.2102997) q[3];
sx q[3];
rz(2.4909693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.66854746) q[2];
sx q[2];
rz(-1.0205525) q[2];
sx q[2];
rz(-2.891053) q[2];
rz(-2.1446832) q[3];
sx q[3];
rz(-1.1397811) q[3];
sx q[3];
rz(0.38058773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32960358) q[0];
sx q[0];
rz(-2.6530837) q[0];
sx q[0];
rz(1.7937775) q[0];
rz(-0.40428058) q[1];
sx q[1];
rz(-0.75899044) q[1];
sx q[1];
rz(-2.0124729) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30325748) q[0];
sx q[0];
rz(-1.1764488) q[0];
sx q[0];
rz(-2.3189937) q[0];
rz(-pi) q[1];
rz(-2.2363844) q[2];
sx q[2];
rz(-1.7844177) q[2];
sx q[2];
rz(-2.7366432) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.68539219) q[1];
sx q[1];
rz(-2.2397352) q[1];
sx q[1];
rz(-0.13369707) q[1];
rz(-pi) q[2];
rz(1.681116) q[3];
sx q[3];
rz(-1.2577783) q[3];
sx q[3];
rz(-0.68389326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3471442) q[2];
sx q[2];
rz(-0.930154) q[2];
sx q[2];
rz(-1.244119) q[2];
rz(2.0813023) q[3];
sx q[3];
rz(-2.5131707) q[3];
sx q[3];
rz(2.8384143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65619549) q[0];
sx q[0];
rz(-3.0363016) q[0];
sx q[0];
rz(-2.7144077) q[0];
rz(-3.1071013) q[1];
sx q[1];
rz(-1.5723615) q[1];
sx q[1];
rz(-3.131386) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1447574) q[0];
sx q[0];
rz(-0.0029276927) q[0];
sx q[0];
rz(2.5064431) q[0];
x q[1];
rz(-0.9385061) q[2];
sx q[2];
rz(-1.7006497) q[2];
sx q[2];
rz(-2.1406895) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2040703) q[1];
sx q[1];
rz(-2.2559705) q[1];
sx q[1];
rz(-3.0584945) q[1];
x q[2];
rz(2.5764546) q[3];
sx q[3];
rz(-1.1332773) q[3];
sx q[3];
rz(1.138162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5222142) q[2];
sx q[2];
rz(-2.7284315) q[2];
sx q[2];
rz(-0.79356066) q[2];
rz(-2.2788952) q[3];
sx q[3];
rz(-1.5669275) q[3];
sx q[3];
rz(2.4376552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4392387) q[0];
sx q[0];
rz(-2.2965501) q[0];
sx q[0];
rz(2.0080361) q[0];
rz(-1.0776862) q[1];
sx q[1];
rz(-0.51529854) q[1];
sx q[1];
rz(2.8394707) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7493499) q[0];
sx q[0];
rz(-1.8729405) q[0];
sx q[0];
rz(-2.7878615) q[0];
rz(-pi) q[1];
rz(2.5100915) q[2];
sx q[2];
rz(-2.4417447) q[2];
sx q[2];
rz(-1.8818784) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4430284) q[1];
sx q[1];
rz(-0.57293597) q[1];
sx q[1];
rz(-2.7354476) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98432912) q[3];
sx q[3];
rz(-0.97198379) q[3];
sx q[3];
rz(-0.58517712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2332396) q[2];
sx q[2];
rz(-0.88963228) q[2];
sx q[2];
rz(1.7944149) q[2];
rz(1.2498648) q[3];
sx q[3];
rz(-1.1976539) q[3];
sx q[3];
rz(-0.10716042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5477448) q[0];
sx q[0];
rz(-2.7084454) q[0];
sx q[0];
rz(-0.10840848) q[0];
rz(-1.0477061) q[1];
sx q[1];
rz(-2.3240418) q[1];
sx q[1];
rz(1.0188867) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86025809) q[0];
sx q[0];
rz(-1.0509914) q[0];
sx q[0];
rz(-2.5297574) q[0];
x q[1];
rz(0.0068130612) q[2];
sx q[2];
rz(-0.76131135) q[2];
sx q[2];
rz(-2.6369264) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0483886) q[1];
sx q[1];
rz(-1.1133514) q[1];
sx q[1];
rz(2.0162575) q[1];
x q[2];
rz(1.4298444) q[3];
sx q[3];
rz(-2.3059855) q[3];
sx q[3];
rz(0.47675374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.63742796) q[2];
sx q[2];
rz(-2.1647858) q[2];
sx q[2];
rz(-2.7039779) q[2];
rz(-2.6484683) q[3];
sx q[3];
rz(-2.2212641) q[3];
sx q[3];
rz(-1.5776618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.86579943) q[0];
sx q[0];
rz(-1.9075305) q[0];
sx q[0];
rz(0.27780521) q[0];
rz(-1.5996784) q[1];
sx q[1];
rz(-2.1733687) q[1];
sx q[1];
rz(0.42617282) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9040462) q[0];
sx q[0];
rz(-1.5874366) q[0];
sx q[0];
rz(0.017701935) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.8502281) q[2];
sx q[2];
rz(-2.2022708) q[2];
sx q[2];
rz(-0.15635083) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3580035) q[1];
sx q[1];
rz(-1.3567909) q[1];
sx q[1];
rz(-3.0215279) q[1];
x q[2];
rz(2.3087193) q[3];
sx q[3];
rz(-2.1578721) q[3];
sx q[3];
rz(-0.47490109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.49491945) q[2];
sx q[2];
rz(-2.9743331) q[2];
sx q[2];
rz(-2.0885928) q[2];
rz(0.35887512) q[3];
sx q[3];
rz(-1.5263298) q[3];
sx q[3];
rz(2.0487962) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36295715) q[0];
sx q[0];
rz(-0.21720049) q[0];
sx q[0];
rz(0.92078513) q[0];
rz(-2.6329363) q[1];
sx q[1];
rz(-2.3743036) q[1];
sx q[1];
rz(2.7210534) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9284436) q[0];
sx q[0];
rz(-2.2290813) q[0];
sx q[0];
rz(0.64853213) q[0];
rz(-1.8968209) q[2];
sx q[2];
rz(-0.65648001) q[2];
sx q[2];
rz(-0.028353779) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1433761) q[1];
sx q[1];
rz(-2.3399379) q[1];
sx q[1];
rz(-2.2364535) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2617662) q[3];
sx q[3];
rz(-1.4925957) q[3];
sx q[3];
rz(1.2134589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.73021182) q[2];
sx q[2];
rz(-2.3513942) q[2];
sx q[2];
rz(-0.31663695) q[2];
rz(0.5591048) q[3];
sx q[3];
rz(-2.5059301) q[3];
sx q[3];
rz(-0.81812286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6560787) q[0];
sx q[0];
rz(-2.06388) q[0];
sx q[0];
rz(-2.0583454) q[0];
rz(0.47646933) q[1];
sx q[1];
rz(-1.0404027) q[1];
sx q[1];
rz(-1.7146005) q[1];
rz(-3.0546679) q[2];
sx q[2];
rz(-2.3200421) q[2];
sx q[2];
rz(2.5210597) q[2];
rz(3.0477897) q[3];
sx q[3];
rz(-1.9524912) q[3];
sx q[3];
rz(-2.5759202) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
