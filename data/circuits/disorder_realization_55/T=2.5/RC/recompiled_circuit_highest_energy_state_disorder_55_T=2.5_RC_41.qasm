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
rz(-1.363938) q[0];
sx q[0];
rz(-2.7231556) q[0];
sx q[0];
rz(1.3016181) q[0];
rz(-2.9786181) q[1];
sx q[1];
rz(-1.6956704) q[1];
sx q[1];
rz(3.1248098) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7201286) q[0];
sx q[0];
rz(-1.6454433) q[0];
sx q[0];
rz(1.9326841) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6616102) q[2];
sx q[2];
rz(-2.9613284) q[2];
sx q[2];
rz(-1.4992876) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.77474817) q[1];
sx q[1];
rz(-1.5758744) q[1];
sx q[1];
rz(-3.1366411) q[1];
x q[2];
rz(0.34256012) q[3];
sx q[3];
rz(-2.9800219) q[3];
sx q[3];
rz(-1.8752961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5590543) q[2];
sx q[2];
rz(-2.4808919) q[2];
sx q[2];
rz(1.5622697) q[2];
rz(-0.2114547) q[3];
sx q[3];
rz(-0.00051694218) q[3];
sx q[3];
rz(2.9805984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.6120537) q[0];
sx q[0];
rz(-2.8654629) q[0];
sx q[0];
rz(1.8147234) q[0];
rz(-0.57948411) q[1];
sx q[1];
rz(-0.0038298413) q[1];
sx q[1];
rz(2.5025867) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2761573) q[0];
sx q[0];
rz(-2.1281181) q[0];
sx q[0];
rz(-2.373567) q[0];
rz(3.0045549) q[2];
sx q[2];
rz(-0.1220905) q[2];
sx q[2];
rz(-0.11954319) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7182448) q[1];
sx q[1];
rz(-3.1210174) q[1];
sx q[1];
rz(-0.59845509) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.78044807) q[3];
sx q[3];
rz(-1.6326346) q[3];
sx q[3];
rz(1.058418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4352033) q[2];
sx q[2];
rz(-3.0053164) q[2];
sx q[2];
rz(-1.603568) q[2];
rz(1.5806574) q[3];
sx q[3];
rz(-3.1272562) q[3];
sx q[3];
rz(0.030979009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.2494217) q[0];
sx q[0];
rz(-2.6264661) q[0];
sx q[0];
rz(-2.7601335) q[0];
rz(-2.4341266) q[1];
sx q[1];
rz(-0.019376945) q[1];
sx q[1];
rz(-1.1245419) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3113925) q[0];
sx q[0];
rz(-1.8233577) q[0];
sx q[0];
rz(0.24858944) q[0];
rz(-pi) q[1];
rz(0.025990268) q[2];
sx q[2];
rz(-1.455869) q[2];
sx q[2];
rz(-2.9851802) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8360236) q[1];
sx q[1];
rz(-1.5587121) q[1];
sx q[1];
rz(-0.062968465) q[1];
rz(2.400983) q[3];
sx q[3];
rz(-0.69878529) q[3];
sx q[3];
rz(1.6388338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6996998) q[2];
sx q[2];
rz(-0.012233891) q[2];
sx q[2];
rz(3.0921248) q[2];
rz(2.5345645) q[3];
sx q[3];
rz(-0.0012461239) q[3];
sx q[3];
rz(-1.2073257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1599051) q[0];
sx q[0];
rz(-0.1683546) q[0];
sx q[0];
rz(3.1244151) q[0];
rz(0.29302868) q[1];
sx q[1];
rz(-2.3511062) q[1];
sx q[1];
rz(1.5944098) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7538504) q[0];
sx q[0];
rz(-1.004129) q[0];
sx q[0];
rz(0.64105861) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3457005) q[2];
sx q[2];
rz(-2.101311) q[2];
sx q[2];
rz(-0.38166416) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.80611011) q[1];
sx q[1];
rz(-1.4467738) q[1];
sx q[1];
rz(-1.6308484) q[1];
rz(2.4956216) q[3];
sx q[3];
rz(-1.5785909) q[3];
sx q[3];
rz(-0.92168671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.25345099) q[2];
sx q[2];
rz(-2.6710822) q[2];
sx q[2];
rz(-0.68224254) q[2];
rz(-3.0864129) q[3];
sx q[3];
rz(-3.1339055) q[3];
sx q[3];
rz(1.3050219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3799915) q[0];
sx q[0];
rz(-0.43552265) q[0];
sx q[0];
rz(-2.7165661) q[0];
rz(1.60166) q[1];
sx q[1];
rz(-0.48307499) q[1];
sx q[1];
rz(2.3262598) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3638301) q[0];
sx q[0];
rz(-1.5095995) q[0];
sx q[0];
rz(-1.5817002) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.710395) q[2];
sx q[2];
rz(-0.023471467) q[2];
sx q[2];
rz(-1.0271629) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7338431) q[1];
sx q[1];
rz(-1.4516648) q[1];
sx q[1];
rz(-1.6581737) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6922705) q[3];
sx q[3];
rz(-1.2290579) q[3];
sx q[3];
rz(-0.69167826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0067979) q[2];
sx q[2];
rz(-0.012601348) q[2];
sx q[2];
rz(1.6689782) q[2];
rz(-0.93004477) q[3];
sx q[3];
rz(-3.127122) q[3];
sx q[3];
rz(0.84021935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8188266) q[0];
sx q[0];
rz(-0.023094026) q[0];
sx q[0];
rz(1.7148788) q[0];
rz(-2.4115883) q[1];
sx q[1];
rz(-0.58861029) q[1];
sx q[1];
rz(1.1013365) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8555785) q[0];
sx q[0];
rz(-1.8402303) q[0];
sx q[0];
rz(2.2982909) q[0];
x q[1];
rz(-0.16898245) q[2];
sx q[2];
rz(-1.797953) q[2];
sx q[2];
rz(2.3262466) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.37820617) q[1];
sx q[1];
rz(-0.15608938) q[1];
sx q[1];
rz(0.915145) q[1];
rz(-pi) q[2];
rz(3.1093842) q[3];
sx q[3];
rz(-1.695249) q[3];
sx q[3];
rz(-2.2197753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.71658984) q[2];
sx q[2];
rz(-3.0807107) q[2];
sx q[2];
rz(-1.3072183) q[2];
rz(0.37846765) q[3];
sx q[3];
rz(-3.1186447) q[3];
sx q[3];
rz(0.65346658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3017479) q[0];
sx q[0];
rz(-1.2675588) q[0];
sx q[0];
rz(0.92754716) q[0];
rz(1.357366) q[1];
sx q[1];
rz(-2.3097242) q[1];
sx q[1];
rz(1.5313139) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6730225) q[0];
sx q[0];
rz(-0.044402145) q[0];
sx q[0];
rz(-1.0442249) q[0];
x q[1];
rz(2.6035735) q[2];
sx q[2];
rz(-2.6876861) q[2];
sx q[2];
rz(1.5858142) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.011259638) q[1];
sx q[1];
rz(-1.5681055) q[1];
sx q[1];
rz(1.4417159) q[1];
x q[2];
rz(-1.0449991) q[3];
sx q[3];
rz(-2.5657095) q[3];
sx q[3];
rz(0.27920846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.36074582) q[2];
sx q[2];
rz(-3.1372034) q[2];
sx q[2];
rz(1.8878262) q[2];
rz(2.4474261) q[3];
sx q[3];
rz(-0.73918754) q[3];
sx q[3];
rz(2.9085801) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53741443) q[0];
sx q[0];
rz(-2.1359213) q[0];
sx q[0];
rz(2.1042714) q[0];
rz(1.6088156) q[1];
sx q[1];
rz(-0.2205801) q[1];
sx q[1];
rz(1.467009) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1572185) q[0];
sx q[0];
rz(-1.3757924) q[0];
sx q[0];
rz(0.068679811) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.59378) q[2];
sx q[2];
rz(-1.8488374) q[2];
sx q[2];
rz(2.8756623) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8716177) q[1];
sx q[1];
rz(-1.570357) q[1];
sx q[1];
rz(1.5719218) q[1];
rz(-pi) q[2];
rz(2.9503533) q[3];
sx q[3];
rz(-2.6789224) q[3];
sx q[3];
rz(1.3706051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.011270114) q[2];
sx q[2];
rz(-2.932817) q[2];
sx q[2];
rz(0.043896349) q[2];
rz(0.52055001) q[3];
sx q[3];
rz(-0.0046516727) q[3];
sx q[3];
rz(-1.6827778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7875824) q[0];
sx q[0];
rz(-3.1391322) q[0];
sx q[0];
rz(-1.3186697) q[0];
rz(1.4175381) q[1];
sx q[1];
rz(-2.8520165) q[1];
sx q[1];
rz(-1.5444548) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.441576) q[0];
sx q[0];
rz(-2.0487983) q[0];
sx q[0];
rz(1.4411323) q[0];
rz(-pi) q[1];
rz(-0.6575281) q[2];
sx q[2];
rz(-1.3597466) q[2];
sx q[2];
rz(-1.576265) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.456157) q[1];
sx q[1];
rz(-0.36000571) q[1];
sx q[1];
rz(2.092157) q[1];
x q[2];
rz(0.61267743) q[3];
sx q[3];
rz(-2.3542488) q[3];
sx q[3];
rz(1.2932216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.80457193) q[2];
sx q[2];
rz(-1.2995517) q[2];
sx q[2];
rz(0.19680944) q[2];
rz(1.9491516) q[3];
sx q[3];
rz(-0.20659031) q[3];
sx q[3];
rz(2.9406252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
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
rz(2.8404959) q[0];
sx q[0];
rz(-1.7844642) q[0];
sx q[0];
rz(-1.1916196) q[0];
rz(-1.5246897) q[1];
sx q[1];
rz(-0.646851) q[1];
sx q[1];
rz(1.5651388) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8904306) q[0];
sx q[0];
rz(-1.403247) q[0];
sx q[0];
rz(-0.29775374) q[0];
x q[1];
rz(-1.1249816) q[2];
sx q[2];
rz(-2.6458394) q[2];
sx q[2];
rz(-3.1204566) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8305739) q[1];
sx q[1];
rz(-1.569848) q[1];
sx q[1];
rz(1.5702973) q[1];
rz(-pi) q[2];
rz(-2.2433167) q[3];
sx q[3];
rz(-0.18096033) q[3];
sx q[3];
rz(-0.74599904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.18371753) q[2];
sx q[2];
rz(-0.58791939) q[2];
sx q[2];
rz(-1.4584165) q[2];
rz(0.030473907) q[3];
sx q[3];
rz(-3.1320429) q[3];
sx q[3];
rz(-2.9392346) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4771582) q[0];
sx q[0];
rz(-1.8071334) q[0];
sx q[0];
rz(-1.4596756) q[0];
rz(1.5674113) q[1];
sx q[1];
rz(-1.3290783) q[1];
sx q[1];
rz(-3.0507416) q[1];
rz(1.6363999) q[2];
sx q[2];
rz(-0.075724307) q[2];
sx q[2];
rz(-2.9157467) q[2];
rz(0.50441691) q[3];
sx q[3];
rz(-0.88263369) q[3];
sx q[3];
rz(0.42019444) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
