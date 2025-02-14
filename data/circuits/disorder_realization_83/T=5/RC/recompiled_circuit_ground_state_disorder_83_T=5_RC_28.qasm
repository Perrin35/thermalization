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
rz(-1.2256844) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8526575) q[0];
sx q[0];
rz(-0.20770099) q[0];
sx q[0];
rz(-0.12909478) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4205409) q[2];
sx q[2];
rz(-1.5051402) q[2];
sx q[2];
rz(-0.70553095) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1154626) q[1];
sx q[1];
rz(-1.5688174) q[1];
sx q[1];
rz(-1.2443365) q[1];
rz(2.1796661) q[3];
sx q[3];
rz(-1.0666218) q[3];
sx q[3];
rz(-3.1162639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.77039069) q[2];
sx q[2];
rz(-2.0870225) q[2];
sx q[2];
rz(-2.2626109) q[2];
rz(-3.0805568) q[3];
sx q[3];
rz(-1.6998484) q[3];
sx q[3];
rz(2.216831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.9691129) q[0];
sx q[0];
rz(-2.0308487) q[0];
sx q[0];
rz(1.2193532) q[0];
rz(-2.1439233) q[1];
sx q[1];
rz(-1.6976633) q[1];
sx q[1];
rz(-2.3692621) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15914842) q[0];
sx q[0];
rz(-1.8799025) q[0];
sx q[0];
rz(-2.3594666) q[0];
rz(-0.017802547) q[2];
sx q[2];
rz(-1.9034198) q[2];
sx q[2];
rz(-2.5359096) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.46281262) q[1];
sx q[1];
rz(-1.2615934) q[1];
sx q[1];
rz(-1.5320381) q[1];
rz(-pi) q[2];
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
rz(2.4274009) q[2];
sx q[2];
rz(-1.6855468) q[2];
sx q[2];
rz(-1.6198772) q[2];
rz(1.6649668) q[3];
sx q[3];
rz(-1.0790389) q[3];
sx q[3];
rz(0.21292201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.49887) q[0];
sx q[0];
rz(-1.1551789) q[0];
sx q[0];
rz(0.89514071) q[0];
rz(2.8816032) q[1];
sx q[1];
rz(-1.7951868) q[1];
sx q[1];
rz(-2.3666429) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74028984) q[0];
sx q[0];
rz(-2.4537114) q[0];
sx q[0];
rz(-2.5581261) q[0];
rz(-1.644448) q[2];
sx q[2];
rz(-1.143874) q[2];
sx q[2];
rz(2.1296453) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5700823) q[1];
sx q[1];
rz(-0.76362839) q[1];
sx q[1];
rz(0.81551733) q[1];
rz(0.633274) q[3];
sx q[3];
rz(-0.97304854) q[3];
sx q[3];
rz(0.55766314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.2669175) q[2];
sx q[2];
rz(-0.83628925) q[2];
sx q[2];
rz(2.3036892) q[2];
rz(-0.72495929) q[3];
sx q[3];
rz(-0.91491142) q[3];
sx q[3];
rz(-0.27423492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2570141) q[0];
sx q[0];
rz(-3.0138636) q[0];
sx q[0];
rz(1.1626441) q[0];
rz(1.1391901) q[1];
sx q[1];
rz(-1.2290686) q[1];
sx q[1];
rz(2.7584279) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1568147) q[0];
sx q[0];
rz(-1.9833282) q[0];
sx q[0];
rz(1.5810313) q[0];
rz(-0.25077925) q[2];
sx q[2];
rz(-2.2927612) q[2];
sx q[2];
rz(2.0112558) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.97505442) q[1];
sx q[1];
rz(-0.34874094) q[1];
sx q[1];
rz(3.1031195) q[1];
rz(-pi) q[2];
rz(-0.64997079) q[3];
sx q[3];
rz(-1.293217) q[3];
sx q[3];
rz(1.7855207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7202516) q[2];
sx q[2];
rz(-0.64695078) q[2];
sx q[2];
rz(1.7163537) q[2];
rz(-1.1228115) q[3];
sx q[3];
rz(-2.4425826) q[3];
sx q[3];
rz(-0.56306806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7177893) q[0];
sx q[0];
rz(-0.20861067) q[0];
sx q[0];
rz(-0.31722379) q[0];
rz(-0.18103655) q[1];
sx q[1];
rz(-1.2174226) q[1];
sx q[1];
rz(-0.6212298) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0962445) q[0];
sx q[0];
rz(-1.3577908) q[0];
sx q[0];
rz(0.81672106) q[0];
rz(-1.6124837) q[2];
sx q[2];
rz(-2.0663849) q[2];
sx q[2];
rz(-2.1245196) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.23074977) q[1];
sx q[1];
rz(-2.158124) q[1];
sx q[1];
rz(-3.0452646) q[1];
x q[2];
rz(2.0192573) q[3];
sx q[3];
rz(-2.3178007) q[3];
sx q[3];
rz(-0.5169249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8165596) q[2];
sx q[2];
rz(-2.4123522) q[2];
sx q[2];
rz(-2.0751591) q[2];
rz(-2.1224497) q[3];
sx q[3];
rz(-0.72503966) q[3];
sx q[3];
rz(-1.9371921) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2524734) q[0];
sx q[0];
rz(-0.42943615) q[0];
sx q[0];
rz(0.47527894) q[0];
rz(-2.1976082) q[1];
sx q[1];
rz(-1.9119268) q[1];
sx q[1];
rz(2.7164187) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19744148) q[0];
sx q[0];
rz(-1.3548491) q[0];
sx q[0];
rz(2.7490389) q[0];
rz(0.026907909) q[2];
sx q[2];
rz(-0.49420658) q[2];
sx q[2];
rz(-1.2721407) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2819967) q[1];
sx q[1];
rz(-2.5841671) q[1];
sx q[1];
rz(1.8962527) q[1];
rz(-pi) q[2];
rz(1.8993072) q[3];
sx q[3];
rz(-1.2385912) q[3];
sx q[3];
rz(-1.1508416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7298594) q[2];
sx q[2];
rz(-1.4889762) q[2];
sx q[2];
rz(-2.497351) q[2];
rz(1.980137) q[3];
sx q[3];
rz(-0.96122733) q[3];
sx q[3];
rz(-0.78288356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2162868) q[0];
sx q[0];
rz(-1.1290978) q[0];
sx q[0];
rz(-2.570545) q[0];
rz(1.1043999) q[1];
sx q[1];
rz(-2.0490502) q[1];
sx q[1];
rz(0.33448321) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12619859) q[0];
sx q[0];
rz(-1.6078464) q[0];
sx q[0];
rz(-1.1657715) q[0];
rz(-1.1715631) q[2];
sx q[2];
rz(-0.81011558) q[2];
sx q[2];
rz(-0.1619815) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.7664288) q[1];
sx q[1];
rz(-0.42459449) q[1];
sx q[1];
rz(0.91318513) q[1];
x q[2];
rz(1.848382) q[3];
sx q[3];
rz(-0.93507877) q[3];
sx q[3];
rz(1.1073974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0337446) q[2];
sx q[2];
rz(-2.0287543) q[2];
sx q[2];
rz(-3.0372078) q[2];
rz(2.8153822) q[3];
sx q[3];
rz(-2.2741337) q[3];
sx q[3];
rz(-0.67702684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0584843) q[0];
sx q[0];
rz(-1.2106189) q[0];
sx q[0];
rz(0.50115681) q[0];
rz(-2.0106563) q[1];
sx q[1];
rz(-0.94291818) q[1];
sx q[1];
rz(1.8556192) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9317497) q[0];
sx q[0];
rz(-2.6476958) q[0];
sx q[0];
rz(-1.2487869) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8752918) q[2];
sx q[2];
rz(-2.2448953) q[2];
sx q[2];
rz(0.69597352) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6567723) q[1];
sx q[1];
rz(-1.3233224) q[1];
sx q[1];
rz(-1.7287068) q[1];
rz(-pi) q[2];
rz(1.0619265) q[3];
sx q[3];
rz(-1.1677907) q[3];
sx q[3];
rz(0.4385192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8274294) q[2];
sx q[2];
rz(-2.8848727) q[2];
sx q[2];
rz(-2.2641505) q[2];
rz(2.1432803) q[3];
sx q[3];
rz(-2.5613997) q[3];
sx q[3];
rz(1.4210526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5494635) q[0];
sx q[0];
rz(-1.2315467) q[0];
sx q[0];
rz(2.8842984) q[0];
rz(-0.65722242) q[1];
sx q[1];
rz(-2.6723599) q[1];
sx q[1];
rz(1.8990272) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22221702) q[0];
sx q[0];
rz(-2.0032481) q[0];
sx q[0];
rz(0.78241703) q[0];
x q[1];
rz(-1.4262582) q[2];
sx q[2];
rz(-0.78510127) q[2];
sx q[2];
rz(0.43971616) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4081289) q[1];
sx q[1];
rz(-2.4545547) q[1];
sx q[1];
rz(-0.71367674) q[1];
x q[2];
rz(-1.7553731) q[3];
sx q[3];
rz(-1.2410302) q[3];
sx q[3];
rz(-2.7884401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2520752) q[2];
sx q[2];
rz(-1.4318848) q[2];
sx q[2];
rz(-0.70402181) q[2];
rz(1.4975123) q[3];
sx q[3];
rz(-1.6234532) q[3];
sx q[3];
rz(1.2161072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7261451) q[0];
sx q[0];
rz(-0.69196597) q[0];
sx q[0];
rz(-0.49219254) q[0];
rz(0.65912143) q[1];
sx q[1];
rz(-0.4159795) q[1];
sx q[1];
rz(0.98512828) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8644442) q[0];
sx q[0];
rz(-1.5168958) q[0];
sx q[0];
rz(0.076907579) q[0];
rz(-pi) q[1];
rz(2.7659723) q[2];
sx q[2];
rz(-1.6723126) q[2];
sx q[2];
rz(0.40292172) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6669504) q[1];
sx q[1];
rz(-1.1792151) q[1];
sx q[1];
rz(-1.8285059) q[1];
rz(-pi) q[2];
rz(-3.1019347) q[3];
sx q[3];
rz(-2.3041906) q[3];
sx q[3];
rz(1.4036251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8681086) q[2];
sx q[2];
rz(-2.2874139) q[2];
sx q[2];
rz(-1.2218529) q[2];
rz(2.6814804) q[3];
sx q[3];
rz(-1.2686138) q[3];
sx q[3];
rz(-2.4251895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7465794) q[0];
sx q[0];
rz(-1.4656675) q[0];
sx q[0];
rz(-1.898484) q[0];
rz(1.6032435) q[1];
sx q[1];
rz(-2.1080882) q[1];
sx q[1];
rz(2.7079667) q[1];
rz(1.4547841) q[2];
sx q[2];
rz(-0.83509904) q[2];
sx q[2];
rz(-0.72088045) q[2];
rz(-3.1258413) q[3];
sx q[3];
rz(-1.2917662) q[3];
sx q[3];
rz(0.36242604) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
