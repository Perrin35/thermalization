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
rz(-0.066601872) q[1];
sx q[1];
rz(5.0636518) q[1];
sx q[1];
rz(10.650462) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9860814) q[0];
sx q[0];
rz(-1.5442463) q[0];
sx q[0];
rz(-2.9355713) q[0];
x q[1];
rz(1.9846693) q[2];
sx q[2];
rz(-2.9777179) q[2];
sx q[2];
rz(-2.6852565) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6027713) q[1];
sx q[1];
rz(-0.32646561) q[1];
sx q[1];
rz(1.5646255) q[1];
rz(-pi) q[2];
rz(-2.1796661) q[3];
sx q[3];
rz(-1.0666218) q[3];
sx q[3];
rz(-0.025328764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.371202) q[2];
sx q[2];
rz(-2.0870225) q[2];
sx q[2];
rz(-0.87898177) q[2];
rz(3.0805568) q[3];
sx q[3];
rz(-1.4417442) q[3];
sx q[3];
rz(2.216831) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1724797) q[0];
sx q[0];
rz(-2.0308487) q[0];
sx q[0];
rz(1.9222395) q[0];
rz(2.1439233) q[1];
sx q[1];
rz(-1.6976633) q[1];
sx q[1];
rz(-0.77233058) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1145086) q[0];
sx q[0];
rz(-0.82875427) q[0];
sx q[0];
rz(-0.42543148) q[0];
x q[1];
rz(-1.2381239) q[2];
sx q[2];
rz(-1.587623) q[2];
sx q[2];
rz(2.1706659) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.58956212) q[1];
sx q[1];
rz(-2.8300474) q[1];
sx q[1];
rz(0.12070848) q[1];
rz(-2.5142383) q[3];
sx q[3];
rz(-1.5511906) q[3];
sx q[3];
rz(-0.62999187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.71419174) q[2];
sx q[2];
rz(-1.6855468) q[2];
sx q[2];
rz(1.5217155) q[2];
rz(1.6649668) q[3];
sx q[3];
rz(-2.0625538) q[3];
sx q[3];
rz(2.9286706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.49887) q[0];
sx q[0];
rz(-1.1551789) q[0];
sx q[0];
rz(2.2464519) q[0];
rz(-0.25998947) q[1];
sx q[1];
rz(-1.7951868) q[1];
sx q[1];
rz(0.77494979) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7827137) q[0];
sx q[0];
rz(-1.2134598) q[0];
sx q[0];
rz(-2.5404447) q[0];
rz(-pi) q[1];
rz(-2.9812212) q[2];
sx q[2];
rz(-2.7087492) q[2];
sx q[2];
rz(-2.3059887) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5436478) q[1];
sx q[1];
rz(-1.0769118) q[1];
sx q[1];
rz(2.1795401) q[1];
rz(-0.633274) q[3];
sx q[3];
rz(-0.97304854) q[3];
sx q[3];
rz(-0.55766314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.2669175) q[2];
sx q[2];
rz(-2.3053034) q[2];
sx q[2];
rz(-2.3036892) q[2];
rz(-2.4166334) q[3];
sx q[3];
rz(-0.91491142) q[3];
sx q[3];
rz(-2.8673577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.88457859) q[0];
sx q[0];
rz(-0.12772904) q[0];
sx q[0];
rz(1.9789486) q[0];
rz(-2.0024025) q[1];
sx q[1];
rz(-1.2290686) q[1];
sx q[1];
rz(2.7584279) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59012198) q[0];
sx q[0];
rz(-1.56142) q[0];
sx q[0];
rz(-0.41255112) q[0];
x q[1];
rz(-1.8454942) q[2];
sx q[2];
rz(-2.3847859) q[2];
sx q[2];
rz(1.5001378) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.63190118) q[1];
sx q[1];
rz(-1.5839403) q[1];
sx q[1];
rz(-2.7930894) q[1];
rz(0.64997079) q[3];
sx q[3];
rz(-1.293217) q[3];
sx q[3];
rz(1.3560719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.42134103) q[2];
sx q[2];
rz(-0.64695078) q[2];
sx q[2];
rz(1.7163537) q[2];
rz(1.1228115) q[3];
sx q[3];
rz(-2.4425826) q[3];
sx q[3];
rz(-2.5785246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4238033) q[0];
sx q[0];
rz(-2.932982) q[0];
sx q[0];
rz(2.8243689) q[0];
rz(2.9605561) q[1];
sx q[1];
rz(-1.2174226) q[1];
sx q[1];
rz(2.5203629) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7214523) q[0];
sx q[0];
rz(-2.3038507) q[0];
sx q[0];
rz(2.8531403) q[0];
rz(-pi) q[1];
x q[1];
rz(0.49595231) q[2];
sx q[2];
rz(-1.5341268) q[2];
sx q[2];
rz(-2.5680361) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9108429) q[1];
sx q[1];
rz(-0.98346868) q[1];
sx q[1];
rz(-0.096328041) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0192573) q[3];
sx q[3];
rz(-2.3178007) q[3];
sx q[3];
rz(-2.6246678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.3250331) q[2];
sx q[2];
rz(-2.4123522) q[2];
sx q[2];
rz(1.0664335) q[2];
rz(2.1224497) q[3];
sx q[3];
rz(-0.72503966) q[3];
sx q[3];
rz(-1.2044005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2524734) q[0];
sx q[0];
rz(-0.42943615) q[0];
sx q[0];
rz(-0.47527894) q[0];
rz(-2.1976082) q[1];
sx q[1];
rz(-1.9119268) q[1];
sx q[1];
rz(-0.42517391) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89579534) q[0];
sx q[0];
rz(-2.6962792) q[0];
sx q[0];
rz(2.6209339) q[0];
rz(-pi) q[1];
rz(2.6475372) q[2];
sx q[2];
rz(-1.5580342) q[2];
sx q[2];
rz(2.8666265) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9035154) q[1];
sx q[1];
rz(-1.0457731) q[1];
sx q[1];
rz(2.9448541) q[1];
x q[2];
rz(-2.7920753) q[3];
sx q[3];
rz(-1.8807285) q[3];
sx q[3];
rz(-2.8323481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4117333) q[2];
sx q[2];
rz(-1.4889762) q[2];
sx q[2];
rz(-2.497351) q[2];
rz(-1.1614557) q[3];
sx q[3];
rz(-2.1803653) q[3];
sx q[3];
rz(0.78288356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92530584) q[0];
sx q[0];
rz(-1.1290978) q[0];
sx q[0];
rz(2.570545) q[0];
rz(2.0371927) q[1];
sx q[1];
rz(-2.0490502) q[1];
sx q[1];
rz(2.8071094) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0153941) q[0];
sx q[0];
rz(-1.6078464) q[0];
sx q[0];
rz(-1.1657715) q[0];
rz(-pi) q[1];
rz(0.38774298) q[2];
sx q[2];
rz(-0.84008145) q[2];
sx q[2];
rz(0.71114388) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4175791) q[1];
sx q[1];
rz(-1.3162606) q[1];
sx q[1];
rz(-1.9144139) q[1];
rz(-2.7860113) q[3];
sx q[3];
rz(-2.4556841) q[3];
sx q[3];
rz(0.65995989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.107848) q[2];
sx q[2];
rz(-1.1128384) q[2];
sx q[2];
rz(3.0372078) q[2];
rz(-0.32621041) q[3];
sx q[3];
rz(-0.86745894) q[3];
sx q[3];
rz(0.67702684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0831083) q[0];
sx q[0];
rz(-1.9309738) q[0];
sx q[0];
rz(-0.50115681) q[0];
rz(-1.1309364) q[1];
sx q[1];
rz(-0.94291818) q[1];
sx q[1];
rz(1.2859734) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0663529) q[0];
sx q[0];
rz(-1.4202002) q[0];
sx q[0];
rz(-2.0429918) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8752918) q[2];
sx q[2];
rz(-2.2448953) q[2];
sx q[2];
rz(-0.69597352) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6567723) q[1];
sx q[1];
rz(-1.3233224) q[1];
sx q[1];
rz(1.4128859) q[1];
rz(-pi) q[2];
rz(2.0796661) q[3];
sx q[3];
rz(-1.1677907) q[3];
sx q[3];
rz(2.7030735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8274294) q[2];
sx q[2];
rz(-0.25671998) q[2];
sx q[2];
rz(2.2641505) q[2];
rz(-2.1432803) q[3];
sx q[3];
rz(-2.5613997) q[3];
sx q[3];
rz(-1.4210526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5921291) q[0];
sx q[0];
rz(-1.9100459) q[0];
sx q[0];
rz(2.8842984) q[0];
rz(-0.65722242) q[1];
sx q[1];
rz(-2.6723599) q[1];
sx q[1];
rz(1.8990272) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9193756) q[0];
sx q[0];
rz(-2.0032481) q[0];
sx q[0];
rz(0.78241703) q[0];
x q[1];
rz(0.14296774) q[2];
sx q[2];
rz(-0.7960628) q[2];
sx q[2];
rz(-0.23676714) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.73346371) q[1];
sx q[1];
rz(-0.68703795) q[1];
sx q[1];
rz(-0.71367674) q[1];
x q[2];
rz(-0.49220645) q[3];
sx q[3];
rz(-2.7653381) q[3];
sx q[3];
rz(0.87615651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.88951748) q[2];
sx q[2];
rz(-1.7097079) q[2];
sx q[2];
rz(2.4375708) q[2];
rz(-1.6440803) q[3];
sx q[3];
rz(-1.5181395) q[3];
sx q[3];
rz(1.9254855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4154476) q[0];
sx q[0];
rz(-0.69196597) q[0];
sx q[0];
rz(0.49219254) q[0];
rz(-0.65912143) q[1];
sx q[1];
rz(-2.7256131) q[1];
sx q[1];
rz(0.98512828) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2378052) q[0];
sx q[0];
rz(-0.093884543) q[0];
sx q[0];
rz(-0.61221377) q[0];
x q[1];
rz(2.7659723) q[2];
sx q[2];
rz(-1.6723126) q[2];
sx q[2];
rz(-2.7386709) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4746423) q[1];
sx q[1];
rz(-1.9623775) q[1];
sx q[1];
rz(1.3130867) q[1];
rz(-pi) q[2];
rz(-0.83701084) q[3];
sx q[3];
rz(-1.5413377) q[3];
sx q[3];
rz(3.000976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8681086) q[2];
sx q[2];
rz(-2.2874139) q[2];
sx q[2];
rz(1.9197397) q[2];
rz(-0.46011225) q[3];
sx q[3];
rz(-1.2686138) q[3];
sx q[3];
rz(0.71640316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39501326) q[0];
sx q[0];
rz(-1.4656675) q[0];
sx q[0];
rz(-1.898484) q[0];
rz(-1.6032435) q[1];
sx q[1];
rz(-1.0335045) q[1];
sx q[1];
rz(-0.43362591) q[1];
rz(3.0144125) q[2];
sx q[2];
rz(-2.3984999) q[2];
sx q[2];
rz(-0.54894757) q[2];
rz(1.5158761) q[3];
sx q[3];
rz(-0.27946278) q[3];
sx q[3];
rz(-2.8362988) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
