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
rz(2.4169154) q[0];
sx q[0];
rz(-1.3114572) q[0];
sx q[0];
rz(0.80479446) q[0];
rz(0.59504396) q[1];
sx q[1];
rz(-2.9437328) q[1];
sx q[1];
rz(-1.9207538) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22053424) q[0];
sx q[0];
rz(-3.0166114) q[0];
sx q[0];
rz(0.28470953) q[0];
x q[1];
rz(0.27313613) q[2];
sx q[2];
rz(-0.47800666) q[2];
sx q[2];
rz(0.82551685) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.28452595) q[1];
sx q[1];
rz(-1.5265122) q[1];
sx q[1];
rz(2.3978049) q[1];
rz(2.9555129) q[3];
sx q[3];
rz(-2.5002989) q[3];
sx q[3];
rz(-2.5945714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7120984) q[2];
sx q[2];
rz(-0.50769371) q[2];
sx q[2];
rz(-1.8053767) q[2];
rz(-1.7441689) q[3];
sx q[3];
rz(-1.5066527) q[3];
sx q[3];
rz(-1.5558745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3812688) q[0];
sx q[0];
rz(-1.5070494) q[0];
sx q[0];
rz(-0.67833483) q[0];
rz(-0.61596576) q[1];
sx q[1];
rz(-2.1149642) q[1];
sx q[1];
rz(-2.9982627) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5472488) q[0];
sx q[0];
rz(-1.7278202) q[0];
sx q[0];
rz(-1.2444522) q[0];
rz(-3.0264261) q[2];
sx q[2];
rz(-0.85571849) q[2];
sx q[2];
rz(-2.5097756) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2611609) q[1];
sx q[1];
rz(-1.5233694) q[1];
sx q[1];
rz(-2.6359184) q[1];
rz(-1.2202699) q[3];
sx q[3];
rz(-1.1683147) q[3];
sx q[3];
rz(0.19798812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4914638) q[2];
sx q[2];
rz(-0.77892196) q[2];
sx q[2];
rz(1.9219363) q[2];
rz(2.8236112) q[3];
sx q[3];
rz(-2.3173083) q[3];
sx q[3];
rz(0.73941755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7431444) q[0];
sx q[0];
rz(-0.92107451) q[0];
sx q[0];
rz(-2.5019124) q[0];
rz(-1.8094874) q[1];
sx q[1];
rz(-2.2001241) q[1];
sx q[1];
rz(2.926362) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9213604) q[0];
sx q[0];
rz(-1.3573933) q[0];
sx q[0];
rz(-2.6409297) q[0];
rz(-pi) q[1];
rz(-1.2864743) q[2];
sx q[2];
rz(-1.5657565) q[2];
sx q[2];
rz(1.7667631) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.97686587) q[1];
sx q[1];
rz(-1.1754719) q[1];
sx q[1];
rz(3.1302451) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3295457) q[3];
sx q[3];
rz(-2.6729306) q[3];
sx q[3];
rz(1.4392343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8362063) q[2];
sx q[2];
rz(-2.6658194) q[2];
sx q[2];
rz(-0.22030182) q[2];
rz(0.78271714) q[3];
sx q[3];
rz(-1.7734776) q[3];
sx q[3];
rz(0.14557423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.43133217) q[0];
sx q[0];
rz(-2.1556957) q[0];
sx q[0];
rz(-1.8248935) q[0];
rz(-0.45007625) q[1];
sx q[1];
rz(-1.4349667) q[1];
sx q[1];
rz(-3.0302474) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5069208) q[0];
sx q[0];
rz(-1.823018) q[0];
sx q[0];
rz(0.66703779) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.088859437) q[2];
sx q[2];
rz(-1.2496762) q[2];
sx q[2];
rz(-0.027983657) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7584472) q[1];
sx q[1];
rz(-2.1229593) q[1];
sx q[1];
rz(0.32917413) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5573978) q[3];
sx q[3];
rz(-1.146493) q[3];
sx q[3];
rz(0.64740023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4075809) q[2];
sx q[2];
rz(-1.13052) q[2];
sx q[2];
rz(-0.13047516) q[2];
rz(2.981251) q[3];
sx q[3];
rz(-0.66157833) q[3];
sx q[3];
rz(-0.17016889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.570356) q[0];
sx q[0];
rz(-2.9354876) q[0];
sx q[0];
rz(-3.0233132) q[0];
rz(-0.89625278) q[1];
sx q[1];
rz(-1.6629442) q[1];
sx q[1];
rz(-0.55312696) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5251585) q[0];
sx q[0];
rz(-1.2431974) q[0];
sx q[0];
rz(2.6382951) q[0];
x q[1];
rz(-0.73605717) q[2];
sx q[2];
rz(-1.2345353) q[2];
sx q[2];
rz(-2.2081809) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70870465) q[1];
sx q[1];
rz(-0.9917534) q[1];
sx q[1];
rz(-2.9396073) q[1];
x q[2];
rz(1.0712264) q[3];
sx q[3];
rz(-2.4132846) q[3];
sx q[3];
rz(2.0508931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7913738) q[2];
sx q[2];
rz(-0.72177902) q[2];
sx q[2];
rz(2.1283894) q[2];
rz(-2.8241217) q[3];
sx q[3];
rz(-1.6976796) q[3];
sx q[3];
rz(2.1875994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1182227) q[0];
sx q[0];
rz(-0.88339266) q[0];
sx q[0];
rz(2.6364442) q[0];
rz(-2.743424) q[1];
sx q[1];
rz(-1.8759556) q[1];
sx q[1];
rz(1.3291976) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5564541) q[0];
sx q[0];
rz(-0.49665652) q[0];
sx q[0];
rz(-1.9209548) q[0];
rz(-pi) q[1];
rz(2.6077725) q[2];
sx q[2];
rz(-1.8835734) q[2];
sx q[2];
rz(-2.52616) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3214282) q[1];
sx q[1];
rz(-1.0099995) q[1];
sx q[1];
rz(3.012077) q[1];
x q[2];
rz(1.8425214) q[3];
sx q[3];
rz(-1.8193805) q[3];
sx q[3];
rz(-2.9523073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.087611) q[2];
sx q[2];
rz(-1.3302646) q[2];
sx q[2];
rz(2.3114253) q[2];
rz(1.4812482) q[3];
sx q[3];
rz(-2.0989428) q[3];
sx q[3];
rz(-1.164485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8692577) q[0];
sx q[0];
rz(-2.2255958) q[0];
sx q[0];
rz(1.2087615) q[0];
rz(-0.2441949) q[1];
sx q[1];
rz(-1.9701651) q[1];
sx q[1];
rz(0.29921439) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9835165) q[0];
sx q[0];
rz(-1.354748) q[0];
sx q[0];
rz(3.1147577) q[0];
rz(0.69912244) q[2];
sx q[2];
rz(-2.4017122) q[2];
sx q[2];
rz(-2.26671) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.077603979) q[1];
sx q[1];
rz(-2.4875459) q[1];
sx q[1];
rz(-1.4908423) q[1];
rz(1.9817623) q[3];
sx q[3];
rz(-2.2310782) q[3];
sx q[3];
rz(-0.11218794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9485335) q[2];
sx q[2];
rz(-1.4024573) q[2];
sx q[2];
rz(-2.194727) q[2];
rz(-1.3777422) q[3];
sx q[3];
rz(-0.54072127) q[3];
sx q[3];
rz(0.64259678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68719012) q[0];
sx q[0];
rz(-2.7070422) q[0];
sx q[0];
rz(0.53892556) q[0];
rz(-0.17351213) q[1];
sx q[1];
rz(-2.3850846) q[1];
sx q[1];
rz(-0.38421806) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.095188524) q[0];
sx q[0];
rz(-2.6676819) q[0];
sx q[0];
rz(-0.35081151) q[0];
x q[1];
rz(-0.60224763) q[2];
sx q[2];
rz(-1.411323) q[2];
sx q[2];
rz(2.1543456) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3620262) q[1];
sx q[1];
rz(-3.0770731) q[1];
sx q[1];
rz(-1.1532591) q[1];
x q[2];
rz(2.124765) q[3];
sx q[3];
rz(-1.5477383) q[3];
sx q[3];
rz(-0.35820828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.82181278) q[2];
sx q[2];
rz(-2.3723448) q[2];
sx q[2];
rz(0.5963076) q[2];
rz(2.1042306) q[3];
sx q[3];
rz(-0.77330247) q[3];
sx q[3];
rz(-1.1843225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7647112) q[0];
sx q[0];
rz(-0.75279555) q[0];
sx q[0];
rz(2.9994614) q[0];
rz(2.0243952) q[1];
sx q[1];
rz(-1.6551207) q[1];
sx q[1];
rz(1.2275009) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1466897) q[0];
sx q[0];
rz(-1.0953383) q[0];
sx q[0];
rz(-1.2593996) q[0];
x q[1];
rz(1.7083211) q[2];
sx q[2];
rz(-1.5548717) q[2];
sx q[2];
rz(-2.8610196) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0490659) q[1];
sx q[1];
rz(-2.5124308) q[1];
sx q[1];
rz(1.4607304) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7015721) q[3];
sx q[3];
rz(-2.2624863) q[3];
sx q[3];
rz(-1.3662076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0907937) q[2];
sx q[2];
rz(-0.24238452) q[2];
sx q[2];
rz(0.79831115) q[2];
rz(1.5797041) q[3];
sx q[3];
rz(-1.3825994) q[3];
sx q[3];
rz(-0.17670512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.056219) q[0];
sx q[0];
rz(-0.64301411) q[0];
sx q[0];
rz(-0.0052848919) q[0];
rz(3.058744) q[1];
sx q[1];
rz(-2.0382035) q[1];
sx q[1];
rz(-2.8740035) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9339211) q[0];
sx q[0];
rz(-0.90822847) q[0];
sx q[0];
rz(1.1712058) q[0];
x q[1];
rz(0.93689367) q[2];
sx q[2];
rz(-2.6724901) q[2];
sx q[2];
rz(-1.7696976) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.96331099) q[1];
sx q[1];
rz(-1.886459) q[1];
sx q[1];
rz(-2.9166978) q[1];
x q[2];
rz(-0.83290108) q[3];
sx q[3];
rz(-0.65147984) q[3];
sx q[3];
rz(-0.92361084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.573632) q[2];
sx q[2];
rz(-0.48144123) q[2];
sx q[2];
rz(1.9285704) q[2];
rz(1.8968449) q[3];
sx q[3];
rz(-2.6601578) q[3];
sx q[3];
rz(1.1904233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5382814) q[0];
sx q[0];
rz(-0.9700226) q[0];
sx q[0];
rz(0.097401311) q[0];
rz(-3.1311323) q[1];
sx q[1];
rz(-2.2757826) q[1];
sx q[1];
rz(-0.59715685) q[1];
rz(0.53562989) q[2];
sx q[2];
rz(-1.4944854) q[2];
sx q[2];
rz(0.22345389) q[2];
rz(1.8298772) q[3];
sx q[3];
rz(-0.6049317) q[3];
sx q[3];
rz(-2.3243454) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
