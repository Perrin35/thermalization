OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9187501) q[0];
sx q[0];
rz(-1.1981244) q[0];
sx q[0];
rz(-2.591748) q[0];
rz(3.0749908) q[1];
sx q[1];
rz(-1.9220592) q[1];
sx q[1];
rz(1.9159082) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15551122) q[0];
sx q[0];
rz(-1.5442463) q[0];
sx q[0];
rz(2.9355713) q[0];
x q[1];
rz(-1.9846693) q[2];
sx q[2];
rz(-2.9777179) q[2];
sx q[2];
rz(-0.45633612) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6027713) q[1];
sx q[1];
rz(-0.32646561) q[1];
sx q[1];
rz(-1.5769671) q[1];
x q[2];
rz(0.80337779) q[3];
sx q[3];
rz(-2.3720494) q[3];
sx q[3];
rz(-0.93955296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77039069) q[2];
sx q[2];
rz(-1.0545701) q[2];
sx q[2];
rz(2.2626109) q[2];
rz(0.061035872) q[3];
sx q[3];
rz(-1.4417442) q[3];
sx q[3];
rz(-2.216831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9691129) q[0];
sx q[0];
rz(-2.0308487) q[0];
sx q[0];
rz(1.2193532) q[0];
rz(-0.99766937) q[1];
sx q[1];
rz(-1.6976633) q[1];
sx q[1];
rz(2.3692621) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9824442) q[0];
sx q[0];
rz(-1.8799025) q[0];
sx q[0];
rz(-0.7821261) q[0];
x q[1];
rz(1.6222811) q[2];
sx q[2];
rz(-2.8085108) q[2];
sx q[2];
rz(-2.5903828) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0454084) q[1];
sx q[1];
rz(-1.6077157) q[1];
sx q[1];
rz(0.30942076) q[1];
x q[2];
rz(-3.1082013) q[3];
sx q[3];
rz(-0.62761939) q[3];
sx q[3];
rz(-2.1737568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4274009) q[2];
sx q[2];
rz(-1.6855468) q[2];
sx q[2];
rz(1.6198772) q[2];
rz(1.6649668) q[3];
sx q[3];
rz(-1.0790389) q[3];
sx q[3];
rz(-2.9286706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.49887) q[0];
sx q[0];
rz(-1.9864137) q[0];
sx q[0];
rz(0.89514071) q[0];
rz(0.25998947) q[1];
sx q[1];
rz(-1.3464059) q[1];
sx q[1];
rz(0.77494979) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4473559) q[0];
sx q[0];
rz(-2.1292344) q[0];
sx q[0];
rz(-1.9959227) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4971447) q[2];
sx q[2];
rz(-1.143874) q[2];
sx q[2];
rz(-2.1296453) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4878825) q[1];
sx q[1];
rz(-1.0431494) q[1];
sx q[1];
rz(0.58077537) q[1];
rz(2.2863359) q[3];
sx q[3];
rz(-2.3000613) q[3];
sx q[3];
rz(-0.35910142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8746752) q[2];
sx q[2];
rz(-2.3053034) q[2];
sx q[2];
rz(-0.8379035) q[2];
rz(-0.72495929) q[3];
sx q[3];
rz(-2.2266812) q[3];
sx q[3];
rz(0.27423492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88457859) q[0];
sx q[0];
rz(-3.0138636) q[0];
sx q[0];
rz(1.9789486) q[0];
rz(-2.0024025) q[1];
sx q[1];
rz(-1.2290686) q[1];
sx q[1];
rz(2.7584279) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1823381) q[0];
sx q[0];
rz(-2.7289411) q[0];
sx q[0];
rz(-3.1182108) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8454942) q[2];
sx q[2];
rz(-2.3847859) q[2];
sx q[2];
rz(-1.5001378) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.93411968) q[1];
sx q[1];
rz(-1.9192682) q[1];
sx q[1];
rz(-1.5568118) q[1];
rz(-pi) q[2];
rz(-1.2270893) q[3];
sx q[3];
rz(-0.94961221) q[3];
sx q[3];
rz(-2.7215001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7202516) q[2];
sx q[2];
rz(-2.4946419) q[2];
sx q[2];
rz(-1.7163537) q[2];
rz(-1.1228115) q[3];
sx q[3];
rz(-0.69901005) q[3];
sx q[3];
rz(0.56306806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7177893) q[0];
sx q[0];
rz(-2.932982) q[0];
sx q[0];
rz(-2.8243689) q[0];
rz(0.18103655) q[1];
sx q[1];
rz(-1.92417) q[1];
sx q[1];
rz(-0.6212298) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.045348195) q[0];
sx q[0];
rz(-1.7838019) q[0];
sx q[0];
rz(-2.3248716) q[0];
rz(-2.6456403) q[2];
sx q[2];
rz(-1.5341268) q[2];
sx q[2];
rz(0.57355659) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8550398) q[1];
sx q[1];
rz(-1.4906487) q[1];
sx q[1];
rz(-2.1602693) q[1];
x q[2];
rz(-2.703691) q[3];
sx q[3];
rz(-0.84841484) q[3];
sx q[3];
rz(-3.0423328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8165596) q[2];
sx q[2];
rz(-2.4123522) q[2];
sx q[2];
rz(2.0751591) q[2];
rz(-2.1224497) q[3];
sx q[3];
rz(-0.72503966) q[3];
sx q[3];
rz(-1.9371921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2524734) q[0];
sx q[0];
rz(-0.42943615) q[0];
sx q[0];
rz(2.6663137) q[0];
rz(2.1976082) q[1];
sx q[1];
rz(-1.2296659) q[1];
sx q[1];
rz(2.7164187) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6797513) q[0];
sx q[0];
rz(-1.1878403) q[0];
sx q[0];
rz(-1.8039077) q[0];
rz(1.5852916) q[2];
sx q[2];
rz(-1.0767848) q[2];
sx q[2];
rz(-1.838889) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85959593) q[1];
sx q[1];
rz(-0.55742555) q[1];
sx q[1];
rz(1.8962527) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.389669) q[3];
sx q[3];
rz(-0.46287352) q[3];
sx q[3];
rz(1.9583256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7298594) q[2];
sx q[2];
rz(-1.6526165) q[2];
sx q[2];
rz(0.64424166) q[2];
rz(1.1614557) q[3];
sx q[3];
rz(-0.96122733) q[3];
sx q[3];
rz(-2.3587091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2162868) q[0];
sx q[0];
rz(-1.1290978) q[0];
sx q[0];
rz(-0.5710477) q[0];
rz(-2.0371927) q[1];
sx q[1];
rz(-1.0925424) q[1];
sx q[1];
rz(-0.33448321) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5307823) q[0];
sx q[0];
rz(-2.73497) q[0];
sx q[0];
rz(-1.4770035) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3399722) q[2];
sx q[2];
rz(-1.2853664) q[2];
sx q[2];
rz(-1.1257671) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7240136) q[1];
sx q[1];
rz(-1.8253321) q[1];
sx q[1];
rz(1.2271787) q[1];
x q[2];
rz(1.2932106) q[3];
sx q[3];
rz(-0.93507877) q[3];
sx q[3];
rz(-1.1073974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0337446) q[2];
sx q[2];
rz(-1.1128384) q[2];
sx q[2];
rz(3.0372078) q[2];
rz(-2.8153822) q[3];
sx q[3];
rz(-0.86745894) q[3];
sx q[3];
rz(-0.67702684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0831083) q[0];
sx q[0];
rz(-1.9309738) q[0];
sx q[0];
rz(-0.50115681) q[0];
rz(-2.0106563) q[1];
sx q[1];
rz(-0.94291818) q[1];
sx q[1];
rz(1.8556192) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0752397) q[0];
sx q[0];
rz(-1.4202002) q[0];
sx q[0];
rz(-1.0986009) q[0];
rz(-pi) q[1];
rz(0.69717631) q[2];
sx q[2];
rz(-1.8072268) q[2];
sx q[2];
rz(-2.4604748) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0946028) q[1];
sx q[1];
rz(-1.7238574) q[1];
sx q[1];
rz(2.8911289) q[1];
rz(-pi) q[2];
rz(0.454161) q[3];
sx q[3];
rz(-1.1060904) q[3];
sx q[3];
rz(-0.91688076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8274294) q[2];
sx q[2];
rz(-2.8848727) q[2];
sx q[2];
rz(2.2641505) q[2];
rz(2.1432803) q[3];
sx q[3];
rz(-2.5613997) q[3];
sx q[3];
rz(1.4210526) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5921291) q[0];
sx q[0];
rz(-1.2315467) q[0];
sx q[0];
rz(-2.8842984) q[0];
rz(-2.4843702) q[1];
sx q[1];
rz(-2.6723599) q[1];
sx q[1];
rz(-1.8990272) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94958828) q[0];
sx q[0];
rz(-2.2704101) q[0];
sx q[0];
rz(-2.5618895) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4262582) q[2];
sx q[2];
rz(-0.78510127) q[2];
sx q[2];
rz(2.7018765) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5661652) q[1];
sx q[1];
rz(-2.0708443) q[1];
sx q[1];
rz(2.0636255) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.49220645) q[3];
sx q[3];
rz(-0.37625458) q[3];
sx q[3];
rz(2.2654361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2520752) q[2];
sx q[2];
rz(-1.4318848) q[2];
sx q[2];
rz(-0.70402181) q[2];
rz(1.6440803) q[3];
sx q[3];
rz(-1.6234532) q[3];
sx q[3];
rz(1.9254855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4154476) q[0];
sx q[0];
rz(-2.4496267) q[0];
sx q[0];
rz(0.49219254) q[0];
rz(0.65912143) q[1];
sx q[1];
rz(-0.4159795) q[1];
sx q[1];
rz(-2.1564644) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2378052) q[0];
sx q[0];
rz(-3.0477081) q[0];
sx q[0];
rz(-2.5293789) q[0];
rz(-pi) q[1];
rz(-0.3756204) q[2];
sx q[2];
rz(-1.4692801) q[2];
sx q[2];
rz(-0.40292172) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6669504) q[1];
sx q[1];
rz(-1.1792151) q[1];
sx q[1];
rz(1.8285059) q[1];
x q[2];
rz(-3.1019347) q[3];
sx q[3];
rz(-0.83740202) q[3];
sx q[3];
rz(-1.4036251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8681086) q[2];
sx q[2];
rz(-0.85417875) q[2];
sx q[2];
rz(-1.2218529) q[2];
rz(-2.6814804) q[3];
sx q[3];
rz(-1.8729788) q[3];
sx q[3];
rz(-2.4251895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.7465794) q[0];
sx q[0];
rz(-1.6759251) q[0];
sx q[0];
rz(1.2431086) q[0];
rz(1.5383491) q[1];
sx q[1];
rz(-1.0335045) q[1];
sx q[1];
rz(-0.43362591) q[1];
rz(-1.4547841) q[2];
sx q[2];
rz(-2.3064936) q[2];
sx q[2];
rz(2.4207122) q[2];
rz(1.6257165) q[3];
sx q[3];
rz(-2.8621299) q[3];
sx q[3];
rz(0.30529387) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
