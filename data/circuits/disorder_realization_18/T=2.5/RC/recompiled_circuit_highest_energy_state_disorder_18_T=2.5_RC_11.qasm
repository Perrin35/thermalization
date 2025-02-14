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
rz(-1.9189605) q[0];
sx q[0];
rz(-0.55813342) q[0];
sx q[0];
rz(-2.4827935) q[0];
rz(0.88762033) q[1];
sx q[1];
rz(2.4467111) q[1];
sx q[1];
rz(12.477787) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6348931) q[0];
sx q[0];
rz(-1.1453712) q[0];
sx q[0];
rz(0.21197196) q[0];
rz(-pi) q[1];
rz(-1.0052571) q[2];
sx q[2];
rz(-1.9449678) q[2];
sx q[2];
rz(-0.94910524) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3693847) q[1];
sx q[1];
rz(-1.7225035) q[1];
sx q[1];
rz(2.1939192) q[1];
x q[2];
rz(-2.1351571) q[3];
sx q[3];
rz(-0.86216037) q[3];
sx q[3];
rz(-0.17911994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6170071) q[2];
sx q[2];
rz(-1.1716537) q[2];
sx q[2];
rz(-0.48801547) q[2];
rz(2.6739142) q[3];
sx q[3];
rz(-2.6165163) q[3];
sx q[3];
rz(0.16744965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92983285) q[0];
sx q[0];
rz(-2.3303895) q[0];
sx q[0];
rz(1.2540586) q[0];
rz(0.60011855) q[1];
sx q[1];
rz(-1.3328726) q[1];
sx q[1];
rz(2.7235203) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.583823) q[0];
sx q[0];
rz(-2.1033786) q[0];
sx q[0];
rz(1.363165) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16440161) q[2];
sx q[2];
rz(-2.4864956) q[2];
sx q[2];
rz(0.10318211) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.046112393) q[1];
sx q[1];
rz(-2.8188092) q[1];
sx q[1];
rz(1.7782636) q[1];
rz(-pi) q[2];
rz(-0.8209867) q[3];
sx q[3];
rz(-2.0625522) q[3];
sx q[3];
rz(-3.0420835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.40689251) q[2];
sx q[2];
rz(-1.1310581) q[2];
sx q[2];
rz(3.0231754) q[2];
rz(-0.40220574) q[3];
sx q[3];
rz(-1.4879358) q[3];
sx q[3];
rz(1.3291298) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6953832) q[0];
sx q[0];
rz(-2.8948247) q[0];
sx q[0];
rz(-1.1545908) q[0];
rz(-2.0788976) q[1];
sx q[1];
rz(-1.0239536) q[1];
sx q[1];
rz(1.1850932) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4615153) q[0];
sx q[0];
rz(-0.74225589) q[0];
sx q[0];
rz(-1.1552325) q[0];
x q[1];
rz(-2.9200771) q[2];
sx q[2];
rz(-2.0755526) q[2];
sx q[2];
rz(-1.5496448) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.94471632) q[1];
sx q[1];
rz(-1.2907952) q[1];
sx q[1];
rz(-1.1796605) q[1];
x q[2];
rz(1.8293158) q[3];
sx q[3];
rz(-1.1791157) q[3];
sx q[3];
rz(1.9323812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9065173) q[2];
sx q[2];
rz(-2.2791028) q[2];
sx q[2];
rz(1.4441351) q[2];
rz(-0.92153543) q[3];
sx q[3];
rz(-2.5369365) q[3];
sx q[3];
rz(-3.0864033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4398572) q[0];
sx q[0];
rz(-0.12049645) q[0];
sx q[0];
rz(0.078027092) q[0];
rz(-2.0241418) q[1];
sx q[1];
rz(-1.8945339) q[1];
sx q[1];
rz(-1.6695581) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3465418) q[0];
sx q[0];
rz(-1.6778089) q[0];
sx q[0];
rz(-2.3404164) q[0];
rz(-2.1536768) q[2];
sx q[2];
rz(-2.0327518) q[2];
sx q[2];
rz(2.5910026) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.3009744) q[1];
sx q[1];
rz(-1.2108016) q[1];
sx q[1];
rz(-0.030351357) q[1];
x q[2];
rz(-1.1999287) q[3];
sx q[3];
rz(-2.2590593) q[3];
sx q[3];
rz(-2.1394068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.83634079) q[2];
sx q[2];
rz(-0.64771104) q[2];
sx q[2];
rz(0.15023896) q[2];
rz(-0.59208208) q[3];
sx q[3];
rz(-1.3034857) q[3];
sx q[3];
rz(-1.1881812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.93064654) q[0];
sx q[0];
rz(-2.7879614) q[0];
sx q[0];
rz(2.3760702) q[0];
rz(0.80134478) q[1];
sx q[1];
rz(-2.0739906) q[1];
sx q[1];
rz(1.9498922) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.649851) q[0];
sx q[0];
rz(-2.5132127) q[0];
sx q[0];
rz(2.6091486) q[0];
rz(-0.52764738) q[2];
sx q[2];
rz(-1.9867267) q[2];
sx q[2];
rz(1.6017583) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5687823) q[1];
sx q[1];
rz(-0.88925075) q[1];
sx q[1];
rz(-3.0589025) q[1];
rz(0.27169835) q[3];
sx q[3];
rz(-1.3110597) q[3];
sx q[3];
rz(1.585497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8754862) q[2];
sx q[2];
rz(-0.85662872) q[2];
sx q[2];
rz(-2.6066656) q[2];
rz(-2.5180425) q[3];
sx q[3];
rz(-2.0303191) q[3];
sx q[3];
rz(2.258544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9542338) q[0];
sx q[0];
rz(-2.7281902) q[0];
sx q[0];
rz(-0.9340539) q[0];
rz(1.2454698) q[1];
sx q[1];
rz(-1.586069) q[1];
sx q[1];
rz(-2.5159871) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68516343) q[0];
sx q[0];
rz(-0.7693839) q[0];
sx q[0];
rz(1.8270232) q[0];
rz(0.90056898) q[2];
sx q[2];
rz(-1.4975542) q[2];
sx q[2];
rz(-1.9894969) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.255068) q[1];
sx q[1];
rz(-2.2400948) q[1];
sx q[1];
rz(-1.943735) q[1];
rz(0.91208338) q[3];
sx q[3];
rz(-0.28277031) q[3];
sx q[3];
rz(-1.2706437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.049456747) q[2];
sx q[2];
rz(-0.9149887) q[2];
sx q[2];
rz(0.12506872) q[2];
rz(-2.4140029) q[3];
sx q[3];
rz(-1.5672507) q[3];
sx q[3];
rz(0.47621134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26559386) q[0];
sx q[0];
rz(-0.47530526) q[0];
sx q[0];
rz(1.884961) q[0];
rz(1.9427293) q[1];
sx q[1];
rz(-1.4224982) q[1];
sx q[1];
rz(1.8667603) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9122353) q[0];
sx q[0];
rz(-1.0692406) q[0];
sx q[0];
rz(-1.0732422) q[0];
rz(-pi) q[1];
rz(0.37097986) q[2];
sx q[2];
rz(-2.0771964) q[2];
sx q[2];
rz(-1.9547403) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3209838) q[1];
sx q[1];
rz(-0.57772355) q[1];
sx q[1];
rz(-2.251365) q[1];
rz(-pi) q[2];
rz(-1.253264) q[3];
sx q[3];
rz(-1.6912563) q[3];
sx q[3];
rz(-2.5250556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4510497) q[2];
sx q[2];
rz(-2.1284911) q[2];
sx q[2];
rz(-0.96735442) q[2];
rz(-1.5585772) q[3];
sx q[3];
rz(-1.6396921) q[3];
sx q[3];
rz(0.30174747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3485182) q[0];
sx q[0];
rz(-2.5418042) q[0];
sx q[0];
rz(3.0317958) q[0];
rz(-1.173136) q[1];
sx q[1];
rz(-0.99183142) q[1];
sx q[1];
rz(-2.0593624) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.968348) q[0];
sx q[0];
rz(-0.81615198) q[0];
sx q[0];
rz(-1.2307172) q[0];
rz(3.033048) q[2];
sx q[2];
rz(-1.0883779) q[2];
sx q[2];
rz(0.71297836) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6746862) q[1];
sx q[1];
rz(-1.162858) q[1];
sx q[1];
rz(0.80927421) q[1];
rz(0.63858648) q[3];
sx q[3];
rz(-1.878384) q[3];
sx q[3];
rz(1.5945376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.75862306) q[2];
sx q[2];
rz(-2.0489645) q[2];
sx q[2];
rz(-0.36337241) q[2];
rz(0.97909561) q[3];
sx q[3];
rz(-0.64371124) q[3];
sx q[3];
rz(0.50160971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8363504) q[0];
sx q[0];
rz(-2.3452121) q[0];
sx q[0];
rz(-2.7889732) q[0];
rz(-1.3615707) q[1];
sx q[1];
rz(-1.9510599) q[1];
sx q[1];
rz(0.1415267) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1175431) q[0];
sx q[0];
rz(-2.96857) q[0];
sx q[0];
rz(-0.19917147) q[0];
rz(-pi) q[1];
rz(0.63371079) q[2];
sx q[2];
rz(-2.4437208) q[2];
sx q[2];
rz(0.77264589) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.334871) q[1];
sx q[1];
rz(-1.1523243) q[1];
sx q[1];
rz(1.970402) q[1];
rz(-pi) q[2];
rz(1.0120506) q[3];
sx q[3];
rz(-2.1390474) q[3];
sx q[3];
rz(-0.38145867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3129468) q[2];
sx q[2];
rz(-1.6673648) q[2];
sx q[2];
rz(-1.0860156) q[2];
rz(-0.18276754) q[3];
sx q[3];
rz(-1.026231) q[3];
sx q[3];
rz(2.1695547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6163841) q[0];
sx q[0];
rz(-1.5807736) q[0];
sx q[0];
rz(0.64424789) q[0];
rz(3.0746025) q[1];
sx q[1];
rz(-1.7552152) q[1];
sx q[1];
rz(-0.11360528) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2841235) q[0];
sx q[0];
rz(-0.22406396) q[0];
sx q[0];
rz(-1.5897008) q[0];
rz(1.468594) q[2];
sx q[2];
rz(-1.488564) q[2];
sx q[2];
rz(-0.036957994) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7305002) q[1];
sx q[1];
rz(-0.20721315) q[1];
sx q[1];
rz(-0.60129182) q[1];
rz(-pi) q[2];
rz(-0.04256256) q[3];
sx q[3];
rz(-1.213338) q[3];
sx q[3];
rz(0.30366158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1429448) q[2];
sx q[2];
rz(-1.8871658) q[2];
sx q[2];
rz(2.3764853) q[2];
rz(0.1782002) q[3];
sx q[3];
rz(-1.5604115) q[3];
sx q[3];
rz(0.83711886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2160303) q[0];
sx q[0];
rz(-2.3860274) q[0];
sx q[0];
rz(-0.26269333) q[0];
rz(-2.8021011) q[1];
sx q[1];
rz(-1.9120293) q[1];
sx q[1];
rz(1.0614352) q[1];
rz(-2.7051666) q[2];
sx q[2];
rz(-0.70311762) q[2];
sx q[2];
rz(1.22618) q[2];
rz(-2.8398561) q[3];
sx q[3];
rz(-1.8724676) q[3];
sx q[3];
rz(-0.68778366) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
