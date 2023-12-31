OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.0032089631) q[0];
sx q[0];
rz(-0.15455833) q[0];
sx q[0];
rz(0.69252339) q[0];
rz(-1.2094296) q[1];
sx q[1];
rz(-1.8930607) q[1];
sx q[1];
rz(-1.7564397) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4614852) q[0];
sx q[0];
rz(-0.84851096) q[0];
sx q[0];
rz(2.0018342) q[0];
rz(-0.42983774) q[2];
sx q[2];
rz(-0.59525437) q[2];
sx q[2];
rz(-2.0855479) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.79917158) q[1];
sx q[1];
rz(-2.85596) q[1];
sx q[1];
rz(-0.61113961) q[1];
rz(-pi) q[2];
rz(0.65269835) q[3];
sx q[3];
rz(-1.1678809) q[3];
sx q[3];
rz(-2.9626915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8866855) q[2];
sx q[2];
rz(-2.343785) q[2];
sx q[2];
rz(0.20516667) q[2];
rz(-0.77130476) q[3];
sx q[3];
rz(-0.78273928) q[3];
sx q[3];
rz(-2.0390959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(2.7339864) q[0];
sx q[0];
rz(-0.74626958) q[0];
sx q[0];
rz(2.6876887) q[0];
rz(1.0247963) q[1];
sx q[1];
rz(-0.4075993) q[1];
sx q[1];
rz(1.9143547) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95943806) q[0];
sx q[0];
rz(-1.6998788) q[0];
sx q[0];
rz(0.076230787) q[0];
rz(2.4279459) q[2];
sx q[2];
rz(-2.4948641) q[2];
sx q[2];
rz(-1.5681859) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.40185803) q[1];
sx q[1];
rz(-1.8622073) q[1];
sx q[1];
rz(1.8406364) q[1];
rz(-pi) q[2];
rz(1.7271872) q[3];
sx q[3];
rz(-1.3723433) q[3];
sx q[3];
rz(-2.1895777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1002905) q[2];
sx q[2];
rz(-1.9561448) q[2];
sx q[2];
rz(-2.5741637) q[2];
rz(-2.7764017) q[3];
sx q[3];
rz(-1.7285715) q[3];
sx q[3];
rz(-2.173483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48297468) q[0];
sx q[0];
rz(-0.56476074) q[0];
sx q[0];
rz(-2.2429402) q[0];
rz(2.1458416) q[1];
sx q[1];
rz(-1.5581222) q[1];
sx q[1];
rz(-0.333289) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62045287) q[0];
sx q[0];
rz(-0.95690173) q[0];
sx q[0];
rz(2.1472473) q[0];
rz(1.0684418) q[2];
sx q[2];
rz(-0.2012673) q[2];
sx q[2];
rz(-1.4301436) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.2644314) q[1];
sx q[1];
rz(-0.96211551) q[1];
sx q[1];
rz(-2.9560637) q[1];
rz(-pi) q[2];
rz(-1.0201449) q[3];
sx q[3];
rz(-1.9198951) q[3];
sx q[3];
rz(0.89494866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4553392) q[2];
sx q[2];
rz(-1.7909966) q[2];
sx q[2];
rz(1.9906445) q[2];
rz(-0.84093705) q[3];
sx q[3];
rz(-2.0261814) q[3];
sx q[3];
rz(1.9870728) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44160098) q[0];
sx q[0];
rz(-1.561152) q[0];
sx q[0];
rz(-2.4568795) q[0];
rz(2.1060064) q[1];
sx q[1];
rz(-2.6338449) q[1];
sx q[1];
rz(1.205014) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7689777) q[0];
sx q[0];
rz(-0.87991558) q[0];
sx q[0];
rz(-3.0983503) q[0];
x q[1];
rz(-0.20385216) q[2];
sx q[2];
rz(-0.9365558) q[2];
sx q[2];
rz(-2.7472592) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.11338621) q[1];
sx q[1];
rz(-0.75062597) q[1];
sx q[1];
rz(1.9764465) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6550001) q[3];
sx q[3];
rz(-2.4349672) q[3];
sx q[3];
rz(-3.1228309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.24923199) q[2];
sx q[2];
rz(-1.426733) q[2];
sx q[2];
rz(-2.7704346) q[2];
rz(1.4012198) q[3];
sx q[3];
rz(-2.4818082) q[3];
sx q[3];
rz(1.1192809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-3.086833) q[0];
sx q[0];
rz(-0.78614569) q[0];
sx q[0];
rz(-3.0084685) q[0];
rz(-0.99331028) q[1];
sx q[1];
rz(-1.3860093) q[1];
sx q[1];
rz(0.55508074) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.426794) q[0];
sx q[0];
rz(-1.5835276) q[0];
sx q[0];
rz(-1.2554332) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.985717) q[2];
sx q[2];
rz(-1.6332111) q[2];
sx q[2];
rz(-1.8780564) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.75025573) q[1];
sx q[1];
rz(-1.5899982) q[1];
sx q[1];
rz(-2.0160497) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6229779) q[3];
sx q[3];
rz(-1.4900041) q[3];
sx q[3];
rz(-2.3063456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.83539) q[2];
sx q[2];
rz(-1.0058879) q[2];
sx q[2];
rz(3.0026657) q[2];
rz(-0.94240087) q[3];
sx q[3];
rz(-1.5706294) q[3];
sx q[3];
rz(0.55148235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54365629) q[0];
sx q[0];
rz(-0.56977001) q[0];
sx q[0];
rz(2.561835) q[0];
rz(-0.12750164) q[1];
sx q[1];
rz(-1.9516877) q[1];
sx q[1];
rz(-1.6019843) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67028763) q[0];
sx q[0];
rz(-1.3677214) q[0];
sx q[0];
rz(2.2905486) q[0];
rz(-0.33467218) q[2];
sx q[2];
rz(-3.0034608) q[2];
sx q[2];
rz(1.7931256) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.74486596) q[1];
sx q[1];
rz(-0.52280451) q[1];
sx q[1];
rz(-0.97739873) q[1];
rz(-1.9220011) q[3];
sx q[3];
rz(-0.71577365) q[3];
sx q[3];
rz(0.90482611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.55398983) q[2];
sx q[2];
rz(-0.25045276) q[2];
sx q[2];
rz(-0.26947752) q[2];
rz(0.23412165) q[3];
sx q[3];
rz(-2.6168489) q[3];
sx q[3];
rz(0.060119303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44678974) q[0];
sx q[0];
rz(-2.5233874) q[0];
sx q[0];
rz(-2.457298) q[0];
rz(-0.11958312) q[1];
sx q[1];
rz(-1.2922829) q[1];
sx q[1];
rz(2.6228242) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25012384) q[0];
sx q[0];
rz(-0.84202535) q[0];
sx q[0];
rz(-0.17983371) q[0];
rz(-pi) q[1];
rz(-1.7905551) q[2];
sx q[2];
rz(-0.81232729) q[2];
sx q[2];
rz(-1.0202368) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.24592933) q[1];
sx q[1];
rz(-2.3024125) q[1];
sx q[1];
rz(0.64114665) q[1];
x q[2];
rz(3.1154237) q[3];
sx q[3];
rz(-1.1702288) q[3];
sx q[3];
rz(-0.27163423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8873022) q[2];
sx q[2];
rz(-1.3672978) q[2];
sx q[2];
rz(2.5781412) q[2];
rz(3.0900132) q[3];
sx q[3];
rz(-1.1392461) q[3];
sx q[3];
rz(-0.95190597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1812487) q[0];
sx q[0];
rz(-0.5287756) q[0];
sx q[0];
rz(1.3990336) q[0];
rz(-2.3545806) q[1];
sx q[1];
rz(-1.1287289) q[1];
sx q[1];
rz(-2.3972437) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3314914) q[0];
sx q[0];
rz(-1.6252675) q[0];
sx q[0];
rz(0.14793747) q[0];
rz(1.0561981) q[2];
sx q[2];
rz(-1.7040164) q[2];
sx q[2];
rz(1.0366057) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.12606049) q[1];
sx q[1];
rz(-2.0704401) q[1];
sx q[1];
rz(1.3152895) q[1];
x q[2];
rz(0.048150496) q[3];
sx q[3];
rz(-0.45504323) q[3];
sx q[3];
rz(-0.20850785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4259592) q[2];
sx q[2];
rz(-1.7461494) q[2];
sx q[2];
rz(0.60950935) q[2];
rz(-0.65731796) q[3];
sx q[3];
rz(-0.64353639) q[3];
sx q[3];
rz(0.26143423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9534) q[0];
sx q[0];
rz(-0.094390079) q[0];
sx q[0];
rz(-1.5040065) q[0];
rz(1.9001182) q[1];
sx q[1];
rz(-1.1499317) q[1];
sx q[1];
rz(2.3666568) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06668815) q[0];
sx q[0];
rz(-0.79830352) q[0];
sx q[0];
rz(1.0134646) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8911792) q[2];
sx q[2];
rz(-0.46386007) q[2];
sx q[2];
rz(-1.8360209) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6248524) q[1];
sx q[1];
rz(-1.7588741) q[1];
sx q[1];
rz(-2.6810758) q[1];
x q[2];
rz(0.998246) q[3];
sx q[3];
rz(-2.8825654) q[3];
sx q[3];
rz(1.5050025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.187414) q[2];
sx q[2];
rz(-2.9118907) q[2];
sx q[2];
rz(0.15787086) q[2];
rz(-1.9291417) q[3];
sx q[3];
rz(-1.0586497) q[3];
sx q[3];
rz(-1.8635748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.071844) q[0];
sx q[0];
rz(-0.97244111) q[0];
sx q[0];
rz(-2.9272595) q[0];
rz(-0.65746039) q[1];
sx q[1];
rz(-2.9174556) q[1];
sx q[1];
rz(-1.0459895) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8911274) q[0];
sx q[0];
rz(-1.2123101) q[0];
sx q[0];
rz(1.6842708) q[0];
x q[1];
rz(-0.59098737) q[2];
sx q[2];
rz(-1.5948442) q[2];
sx q[2];
rz(1.9901333) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9606049) q[1];
sx q[1];
rz(-2.1121896) q[1];
sx q[1];
rz(-2.7156746) q[1];
rz(-pi) q[2];
rz(-0.06185992) q[3];
sx q[3];
rz(-0.40611551) q[3];
sx q[3];
rz(-1.5554242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.41632286) q[2];
sx q[2];
rz(-0.060083397) q[2];
sx q[2];
rz(-2.0521169) q[2];
rz(1.5661092) q[3];
sx q[3];
rz(-1.8982866) q[3];
sx q[3];
rz(-0.65264788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(1.8626704) q[0];
sx q[0];
rz(-0.60386064) q[0];
sx q[0];
rz(0.84558564) q[0];
rz(1.5325585) q[1];
sx q[1];
rz(-1.6747723) q[1];
sx q[1];
rz(2.0369045) q[1];
rz(2.5074742) q[2];
sx q[2];
rz(-2.6478883) q[2];
sx q[2];
rz(-1.4512856) q[2];
rz(-0.10272051) q[3];
sx q[3];
rz(-2.5892047) q[3];
sx q[3];
rz(1.8451167) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
