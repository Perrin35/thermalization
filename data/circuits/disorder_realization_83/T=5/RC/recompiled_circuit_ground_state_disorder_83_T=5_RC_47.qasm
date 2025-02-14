OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(5.0603428) q[0];
sx q[0];
rz(4.3397171) q[0];
sx q[0];
rz(8.8749333) q[0];
rz(3.0749908) q[1];
sx q[1];
rz(-1.9220592) q[1];
sx q[1];
rz(1.9159082) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2889351) q[0];
sx q[0];
rz(-0.20770099) q[0];
sx q[0];
rz(0.12909478) q[0];
x q[1];
rz(1.9846693) q[2];
sx q[2];
rz(-0.16387476) q[2];
sx q[2];
rz(-0.45633612) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1154626) q[1];
sx q[1];
rz(-1.5727753) q[1];
sx q[1];
rz(-1.2443365) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96192653) q[3];
sx q[3];
rz(-2.0749708) q[3];
sx q[3];
rz(-0.025328764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.77039069) q[2];
sx q[2];
rz(-2.0870225) q[2];
sx q[2];
rz(-0.87898177) q[2];
rz(-0.061035872) q[3];
sx q[3];
rz(-1.6998484) q[3];
sx q[3];
rz(0.92476168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9691129) q[0];
sx q[0];
rz(-1.1107439) q[0];
sx q[0];
rz(1.2193532) q[0];
rz(-2.1439233) q[1];
sx q[1];
rz(-1.4439293) q[1];
sx q[1];
rz(-0.77233058) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7051429) q[0];
sx q[0];
rz(-2.3068743) q[0];
sx q[0];
rz(1.993773) q[0];
rz(-pi) q[1];
rz(3.1237901) q[2];
sx q[2];
rz(-1.9034198) q[2];
sx q[2];
rz(-2.5359096) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0454084) q[1];
sx q[1];
rz(-1.533877) q[1];
sx q[1];
rz(0.30942076) q[1];
x q[2];
rz(-0.62735438) q[3];
sx q[3];
rz(-1.590402) q[3];
sx q[3];
rz(2.5116008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4274009) q[2];
sx q[2];
rz(-1.4560459) q[2];
sx q[2];
rz(-1.6198772) q[2];
rz(1.6649668) q[3];
sx q[3];
rz(-1.0790389) q[3];
sx q[3];
rz(-2.9286706) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.49887) q[0];
sx q[0];
rz(-1.9864137) q[0];
sx q[0];
rz(-0.89514071) q[0];
rz(-0.25998947) q[1];
sx q[1];
rz(-1.7951868) q[1];
sx q[1];
rz(0.77494979) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7827137) q[0];
sx q[0];
rz(-1.2134598) q[0];
sx q[0];
rz(0.60114791) q[0];
rz(-pi) q[1];
rz(-1.644448) q[2];
sx q[2];
rz(-1.9977187) q[2];
sx q[2];
rz(-2.1296453) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4878825) q[1];
sx q[1];
rz(-2.0984432) q[1];
sx q[1];
rz(0.58077537) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86943961) q[3];
sx q[3];
rz(-2.0819391) q[3];
sx q[3];
rz(-1.7366586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.2669175) q[2];
sx q[2];
rz(-0.83628925) q[2];
sx q[2];
rz(0.8379035) q[2];
rz(0.72495929) q[3];
sx q[3];
rz(-2.2266812) q[3];
sx q[3];
rz(-0.27423492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88457859) q[0];
sx q[0];
rz(-0.12772904) q[0];
sx q[0];
rz(1.9789486) q[0];
rz(1.1391901) q[1];
sx q[1];
rz(-1.2290686) q[1];
sx q[1];
rz(-0.3831648) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1568147) q[0];
sx q[0];
rz(-1.9833282) q[0];
sx q[0];
rz(1.5605613) q[0];
x q[1];
rz(-2.3085528) q[2];
sx q[2];
rz(-1.383457) q[2];
sx q[2];
rz(0.27275547) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.207473) q[1];
sx q[1];
rz(-1.9192682) q[1];
sx q[1];
rz(-1.5568118) q[1];
x q[2];
rz(1.9145033) q[3];
sx q[3];
rz(-2.1919804) q[3];
sx q[3];
rz(-0.4200926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7202516) q[2];
sx q[2];
rz(-0.64695078) q[2];
sx q[2];
rz(1.7163537) q[2];
rz(-2.0187812) q[3];
sx q[3];
rz(-0.69901005) q[3];
sx q[3];
rz(2.5785246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4238033) q[0];
sx q[0];
rz(-0.20861067) q[0];
sx q[0];
rz(2.8243689) q[0];
rz(0.18103655) q[1];
sx q[1];
rz(-1.2174226) q[1];
sx q[1];
rz(-2.5203629) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3040627) q[0];
sx q[0];
rz(-0.777839) q[0];
sx q[0];
rz(-1.2647948) q[0];
rz(-pi) q[1];
rz(-1.6124837) q[2];
sx q[2];
rz(-1.0752077) q[2];
sx q[2];
rz(2.1245196) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8550398) q[1];
sx q[1];
rz(-1.650944) q[1];
sx q[1];
rz(0.98132332) q[1];
rz(-pi) q[2];
rz(1.1223354) q[3];
sx q[3];
rz(-0.82379195) q[3];
sx q[3];
rz(2.6246678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.3250331) q[2];
sx q[2];
rz(-2.4123522) q[2];
sx q[2];
rz(-2.0751591) q[2];
rz(1.019143) q[3];
sx q[3];
rz(-0.72503966) q[3];
sx q[3];
rz(-1.9371921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8891193) q[0];
sx q[0];
rz(-2.7121565) q[0];
sx q[0];
rz(-0.47527894) q[0];
rz(0.94398445) q[1];
sx q[1];
rz(-1.9119268) q[1];
sx q[1];
rz(-0.42517391) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19744148) q[0];
sx q[0];
rz(-1.7867435) q[0];
sx q[0];
rz(2.7490389) q[0];
rz(-pi) q[1];
x q[1];
rz(1.556301) q[2];
sx q[2];
rz(-2.0648079) q[2];
sx q[2];
rz(-1.838889) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.85959593) q[1];
sx q[1];
rz(-0.55742555) q[1];
sx q[1];
rz(-1.8962527) q[1];
rz(-2.389669) q[3];
sx q[3];
rz(-0.46287352) q[3];
sx q[3];
rz(-1.183267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7298594) q[2];
sx q[2];
rz(-1.6526165) q[2];
sx q[2];
rz(-0.64424166) q[2];
rz(1.1614557) q[3];
sx q[3];
rz(-0.96122733) q[3];
sx q[3];
rz(0.78288356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92530584) q[0];
sx q[0];
rz(-2.0124948) q[0];
sx q[0];
rz(-2.570545) q[0];
rz(-2.0371927) q[1];
sx q[1];
rz(-1.0925424) q[1];
sx q[1];
rz(2.8071094) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0153941) q[0];
sx q[0];
rz(-1.5337463) q[0];
sx q[0];
rz(-1.1657715) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1715631) q[2];
sx q[2];
rz(-0.81011558) q[2];
sx q[2];
rz(0.1619815) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.063364048) q[1];
sx q[1];
rz(-1.9029088) q[1];
sx q[1];
rz(-2.8719877) q[1];
x q[2];
rz(-2.7860113) q[3];
sx q[3];
rz(-0.68590859) q[3];
sx q[3];
rz(2.4816328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.107848) q[2];
sx q[2];
rz(-1.1128384) q[2];
sx q[2];
rz(3.0372078) q[2];
rz(-2.8153822) q[3];
sx q[3];
rz(-2.2741337) q[3];
sx q[3];
rz(0.67702684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0584843) q[0];
sx q[0];
rz(-1.2106189) q[0];
sx q[0];
rz(2.6404358) q[0];
rz(-1.1309364) q[1];
sx q[1];
rz(-0.94291818) q[1];
sx q[1];
rz(-1.8556192) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9317497) q[0];
sx q[0];
rz(-2.6476958) q[0];
sx q[0];
rz(1.8928058) q[0];
x q[1];
rz(1.2663009) q[2];
sx q[2];
rz(-0.89669734) q[2];
sx q[2];
rz(0.69597352) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0803538) q[1];
sx q[1];
rz(-0.29268943) q[1];
sx q[1];
rz(2.5848542) q[1];
rz(-0.85190947) q[3];
sx q[3];
rz(-0.63792568) q[3];
sx q[3];
rz(-1.3965811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3141632) q[2];
sx q[2];
rz(-2.8848727) q[2];
sx q[2];
rz(2.2641505) q[2];
rz(0.99831239) q[3];
sx q[3];
rz(-2.5613997) q[3];
sx q[3];
rz(1.72054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5921291) q[0];
sx q[0];
rz(-1.9100459) q[0];
sx q[0];
rz(-0.25729427) q[0];
rz(0.65722242) q[1];
sx q[1];
rz(-2.6723599) q[1];
sx q[1];
rz(-1.8990272) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94958828) q[0];
sx q[0];
rz(-2.2704101) q[0];
sx q[0];
rz(2.5618895) q[0];
rz(-pi) q[1];
x q[1];
rz(0.79093604) q[2];
sx q[2];
rz(-1.4688014) q[2];
sx q[2];
rz(1.2336709) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8942106) q[1];
sx q[1];
rz(-1.1426434) q[1];
sx q[1];
rz(0.55511554) q[1];
rz(-pi) q[2];
rz(0.49220645) q[3];
sx q[3];
rz(-0.37625458) q[3];
sx q[3];
rz(0.87615651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2520752) q[2];
sx q[2];
rz(-1.4318848) q[2];
sx q[2];
rz(-0.70402181) q[2];
rz(-1.6440803) q[3];
sx q[3];
rz(-1.5181395) q[3];
sx q[3];
rz(1.9254855) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4154476) q[0];
sx q[0];
rz(-2.4496267) q[0];
sx q[0];
rz(0.49219254) q[0];
rz(2.4824712) q[1];
sx q[1];
rz(-2.7256131) q[1];
sx q[1];
rz(0.98512828) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27714849) q[0];
sx q[0];
rz(-1.5168958) q[0];
sx q[0];
rz(3.0646851) q[0];
rz(-2.8707377) q[2];
sx q[2];
rz(-0.388467) q[2];
sx q[2];
rz(2.2252639) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4746423) q[1];
sx q[1];
rz(-1.9623775) q[1];
sx q[1];
rz(-1.3130867) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6147695) q[3];
sx q[3];
rz(-0.73426651) q[3];
sx q[3];
rz(1.6787613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.27348406) q[2];
sx q[2];
rz(-2.2874139) q[2];
sx q[2];
rz(1.9197397) q[2];
rz(-2.6814804) q[3];
sx q[3];
rz(-1.2686138) q[3];
sx q[3];
rz(2.4251895) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7465794) q[0];
sx q[0];
rz(-1.4656675) q[0];
sx q[0];
rz(-1.898484) q[0];
rz(-1.6032435) q[1];
sx q[1];
rz(-1.0335045) q[1];
sx q[1];
rz(-0.43362591) q[1];
rz(1.6868085) q[2];
sx q[2];
rz(-2.3064936) q[2];
sx q[2];
rz(2.4207122) q[2];
rz(-1.6257165) q[3];
sx q[3];
rz(-0.27946278) q[3];
sx q[3];
rz(-2.8362988) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
