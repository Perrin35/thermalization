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
rz(-1.5440829) q[0];
sx q[0];
rz(-1.7680327) q[0];
sx q[0];
rz(-1.0184259) q[0];
rz(-0.51813689) q[1];
sx q[1];
rz(-2.3845446) q[1];
sx q[1];
rz(0.63176027) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1904097) q[0];
sx q[0];
rz(-0.24226878) q[0];
sx q[0];
rz(-2.1184068) q[0];
rz(-0.024876923) q[2];
sx q[2];
rz(-1.8734249) q[2];
sx q[2];
rz(0.51371117) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2871197) q[1];
sx q[1];
rz(-1.2911671) q[1];
sx q[1];
rz(0.92052144) q[1];
rz(1.8381836) q[3];
sx q[3];
rz(-2.8041556) q[3];
sx q[3];
rz(-2.9149559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.53039256) q[2];
sx q[2];
rz(-2.7860614) q[2];
sx q[2];
rz(1.6543039) q[2];
rz(0.90855956) q[3];
sx q[3];
rz(-2.9049951) q[3];
sx q[3];
rz(-2.1366185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0205883) q[0];
sx q[0];
rz(-1.4326743) q[0];
sx q[0];
rz(-2.9793136) q[0];
rz(2.8970215) q[1];
sx q[1];
rz(-1.9385612) q[1];
sx q[1];
rz(2.8499106) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2431716) q[0];
sx q[0];
rz(-0.34472877) q[0];
sx q[0];
rz(0.60144641) q[0];
x q[1];
rz(2.9180718) q[2];
sx q[2];
rz(-1.8475899) q[2];
sx q[2];
rz(-1.8794488) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.32714601) q[1];
sx q[1];
rz(-1.3400048) q[1];
sx q[1];
rz(0.77951851) q[1];
rz(-1.1532182) q[3];
sx q[3];
rz(-1.1912701) q[3];
sx q[3];
rz(2.1714927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.67819277) q[2];
sx q[2];
rz(-2.2726111) q[2];
sx q[2];
rz(1.0758859) q[2];
rz(1.894527) q[3];
sx q[3];
rz(-1.5445292) q[3];
sx q[3];
rz(-1.4472848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4397044) q[0];
sx q[0];
rz(-1.1264369) q[0];
sx q[0];
rz(2.9578748) q[0];
rz(1.6366929) q[1];
sx q[1];
rz(-0.74438649) q[1];
sx q[1];
rz(-0.16477975) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8863109) q[0];
sx q[0];
rz(-1.5802529) q[0];
sx q[0];
rz(-0.9849809) q[0];
rz(0.38305958) q[2];
sx q[2];
rz(-0.7128517) q[2];
sx q[2];
rz(2.0682316) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2307869) q[1];
sx q[1];
rz(-1.9525098) q[1];
sx q[1];
rz(2.4934105) q[1];
x q[2];
rz(2.1042473) q[3];
sx q[3];
rz(-1.8672393) q[3];
sx q[3];
rz(0.25755778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.034721) q[2];
sx q[2];
rz(-1.9811337) q[2];
sx q[2];
rz(1.1064233) q[2];
rz(2.4404081) q[3];
sx q[3];
rz(-1.3283575) q[3];
sx q[3];
rz(-0.4655233) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3727386) q[0];
sx q[0];
rz(-1.1356069) q[0];
sx q[0];
rz(0.94451529) q[0];
rz(1.5185897) q[1];
sx q[1];
rz(-2.1606162) q[1];
sx q[1];
rz(-1.5400344) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2006783) q[0];
sx q[0];
rz(-0.8884065) q[0];
sx q[0];
rz(-2.3290578) q[0];
rz(-pi) q[1];
rz(0.77727274) q[2];
sx q[2];
rz(-0.16675719) q[2];
sx q[2];
rz(0.47698944) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9303927) q[1];
sx q[1];
rz(-0.89387776) q[1];
sx q[1];
rz(-2.6571214) q[1];
rz(-pi) q[2];
rz(-2.3911658) q[3];
sx q[3];
rz(-1.5380757) q[3];
sx q[3];
rz(-2.1145543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1314142) q[2];
sx q[2];
rz(-0.62128908) q[2];
sx q[2];
rz(-0.16769257) q[2];
rz(-0.0532648) q[3];
sx q[3];
rz(-2.0361418) q[3];
sx q[3];
rz(-1.7648599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57529706) q[0];
sx q[0];
rz(-0.10558858) q[0];
sx q[0];
rz(2.8422624) q[0];
rz(-2.7686367) q[1];
sx q[1];
rz(-1.2495709) q[1];
sx q[1];
rz(-1.6711055) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9540967) q[0];
sx q[0];
rz(-2.4135655) q[0];
sx q[0];
rz(2.4690829) q[0];
x q[1];
rz(0.23678933) q[2];
sx q[2];
rz(-2.6663187) q[2];
sx q[2];
rz(1.7101814) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.336323) q[1];
sx q[1];
rz(-1.0674369) q[1];
sx q[1];
rz(0.663228) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3277131) q[3];
sx q[3];
rz(-2.5057) q[3];
sx q[3];
rz(-2.2896374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2028929) q[2];
sx q[2];
rz(-2.4788269) q[2];
sx q[2];
rz(1.1104442) q[2];
rz(-2.2796196) q[3];
sx q[3];
rz(-1.5098666) q[3];
sx q[3];
rz(1.8814258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92457572) q[0];
sx q[0];
rz(-0.68704263) q[0];
sx q[0];
rz(0.4050912) q[0];
rz(-0.336126) q[1];
sx q[1];
rz(-1.7205709) q[1];
sx q[1];
rz(0.95132557) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8444871) q[0];
sx q[0];
rz(-1.9713624) q[0];
sx q[0];
rz(1.8482918) q[0];
x q[1];
rz(-0.23299322) q[2];
sx q[2];
rz(-2.39175) q[2];
sx q[2];
rz(-2.1342333) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8321633) q[1];
sx q[1];
rz(-2.1346666) q[1];
sx q[1];
rz(1.840074) q[1];
x q[2];
rz(-1.7634994) q[3];
sx q[3];
rz(-1.2964222) q[3];
sx q[3];
rz(-1.5128795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.090791) q[2];
sx q[2];
rz(-1.4046706) q[2];
sx q[2];
rz(0.23078272) q[2];
rz(0.18203059) q[3];
sx q[3];
rz(-2.5735276) q[3];
sx q[3];
rz(-2.6128984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7886605) q[0];
sx q[0];
rz(-0.039529888) q[0];
sx q[0];
rz(2.2094862) q[0];
rz(0.63356361) q[1];
sx q[1];
rz(-2.0220058) q[1];
sx q[1];
rz(1.3379785) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0377696) q[0];
sx q[0];
rz(-1.4892206) q[0];
sx q[0];
rz(1.4647116) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4557082) q[2];
sx q[2];
rz(-1.2988) q[2];
sx q[2];
rz(0.65036648) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7576068) q[1];
sx q[1];
rz(-1.9644871) q[1];
sx q[1];
rz(-0.91169731) q[1];
rz(3.0831117) q[3];
sx q[3];
rz(-2.6425458) q[3];
sx q[3];
rz(2.8935695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6771217) q[2];
sx q[2];
rz(-1.2268343) q[2];
sx q[2];
rz(1.9390437) q[2];
rz(-0.36457148) q[3];
sx q[3];
rz(-2.4489844) q[3];
sx q[3];
rz(-2.306849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6637591) q[0];
sx q[0];
rz(-0.57357016) q[0];
sx q[0];
rz(-0.17663503) q[0];
rz(-2.7015576) q[1];
sx q[1];
rz(-1.9562079) q[1];
sx q[1];
rz(0.96493351) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2595183) q[0];
sx q[0];
rz(-1.4384369) q[0];
sx q[0];
rz(2.9455723) q[0];
rz(-pi) q[1];
rz(1.0267443) q[2];
sx q[2];
rz(-2.1026582) q[2];
sx q[2];
rz(0.085467664) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8254556) q[1];
sx q[1];
rz(-1.2082005) q[1];
sx q[1];
rz(2.1522002) q[1];
x q[2];
rz(-3.065751) q[3];
sx q[3];
rz(-1.3170529) q[3];
sx q[3];
rz(2.6395519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9289916) q[2];
sx q[2];
rz(-1.6182199) q[2];
sx q[2];
rz(-3.1316481) q[2];
rz(0.11659226) q[3];
sx q[3];
rz(-2.8023585) q[3];
sx q[3];
rz(-0.54178437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5803439) q[0];
sx q[0];
rz(-1.2224226) q[0];
sx q[0];
rz(0.62193459) q[0];
rz(-1.7033345) q[1];
sx q[1];
rz(-0.57241264) q[1];
sx q[1];
rz(-0.66351801) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2885482) q[0];
sx q[0];
rz(-0.85381928) q[0];
sx q[0];
rz(1.8621481) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.12309317) q[2];
sx q[2];
rz(-1.2893587) q[2];
sx q[2];
rz(1.7116473) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2474644) q[1];
sx q[1];
rz(-0.93527764) q[1];
sx q[1];
rz(2.7429917) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6648977) q[3];
sx q[3];
rz(-1.851436) q[3];
sx q[3];
rz(-0.68583127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5180987) q[2];
sx q[2];
rz(-1.3891209) q[2];
sx q[2];
rz(0.36250472) q[2];
rz(-2.6643961) q[3];
sx q[3];
rz(-1.0537182) q[3];
sx q[3];
rz(-2.0188735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.057673205) q[0];
sx q[0];
rz(-0.73129439) q[0];
sx q[0];
rz(2.2398563) q[0];
rz(-0.67507356) q[1];
sx q[1];
rz(-2.156064) q[1];
sx q[1];
rz(0.14428446) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1146204) q[0];
sx q[0];
rz(-1.9211702) q[0];
sx q[0];
rz(-3.0536985) q[0];
x q[1];
rz(1.6727757) q[2];
sx q[2];
rz(-1.7893409) q[2];
sx q[2];
rz(0.25221014) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.212483) q[1];
sx q[1];
rz(-1.3204201) q[1];
sx q[1];
rz(-0.65156619) q[1];
rz(1.7606401) q[3];
sx q[3];
rz(-2.0060872) q[3];
sx q[3];
rz(-1.3224885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6361864) q[2];
sx q[2];
rz(-1.212684) q[2];
sx q[2];
rz(-0.20067659) q[2];
rz(-2.0319273) q[3];
sx q[3];
rz(-0.4626914) q[3];
sx q[3];
rz(0.5717352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5889482) q[0];
sx q[0];
rz(-2.1727967) q[0];
sx q[0];
rz(-1.1217242) q[0];
rz(2.6523392) q[1];
sx q[1];
rz(-1.6901292) q[1];
sx q[1];
rz(2.1314175) q[1];
rz(0.092838661) q[2];
sx q[2];
rz(-1.1576817) q[2];
sx q[2];
rz(0.67410034) q[2];
rz(1.5420001) q[3];
sx q[3];
rz(-1.1287108) q[3];
sx q[3];
rz(-2.5870833) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
