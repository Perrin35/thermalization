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
rz(0.69005203) q[0];
sx q[0];
rz(-0.85059387) q[0];
sx q[0];
rz(-2.9414862) q[0];
rz(2.9432358) q[1];
sx q[1];
rz(-2.618572) q[1];
sx q[1];
rz(-1.3318292) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49591767) q[0];
sx q[0];
rz(-2.0089275) q[0];
sx q[0];
rz(-2.8854831) q[0];
rz(-pi) q[1];
rz(-0.66352377) q[2];
sx q[2];
rz(-0.71394701) q[2];
sx q[2];
rz(-2.6235142) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7847474) q[1];
sx q[1];
rz(-2.4618399) q[1];
sx q[1];
rz(3.0650782) q[1];
x q[2];
rz(2.0468007) q[3];
sx q[3];
rz(-1.4972357) q[3];
sx q[3];
rz(-1.6122163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0693822) q[2];
sx q[2];
rz(-2.78077) q[2];
sx q[2];
rz(-1.3145831) q[2];
rz(0.9032816) q[3];
sx q[3];
rz(-1.7777781) q[3];
sx q[3];
rz(-0.34590736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2717993) q[0];
sx q[0];
rz(-1.485774) q[0];
sx q[0];
rz(-2.2858009) q[0];
rz(2.3675512) q[1];
sx q[1];
rz(-1.1639405) q[1];
sx q[1];
rz(0.78786293) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5408527) q[0];
sx q[0];
rz(-2.8898281) q[0];
sx q[0];
rz(1.6600079) q[0];
rz(-1.6131975) q[2];
sx q[2];
rz(-1.2268492) q[2];
sx q[2];
rz(2.6478772) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.98709244) q[1];
sx q[1];
rz(-2.5804511) q[1];
sx q[1];
rz(2.2461056) q[1];
rz(-0.26420748) q[3];
sx q[3];
rz(-0.6843066) q[3];
sx q[3];
rz(-0.58240376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0756125) q[2];
sx q[2];
rz(-1.9115261) q[2];
sx q[2];
rz(-0.39609972) q[2];
rz(-1.570545) q[3];
sx q[3];
rz(-1.5038265) q[3];
sx q[3];
rz(1.3554696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9539255) q[0];
sx q[0];
rz(-1.4952156) q[0];
sx q[0];
rz(0.060610108) q[0];
rz(-2.4568779) q[1];
sx q[1];
rz(-1.741332) q[1];
sx q[1];
rz(0.33904591) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44797541) q[0];
sx q[0];
rz(-2.1537499) q[0];
sx q[0];
rz(1.5459803) q[0];
x q[1];
rz(-2.9525063) q[2];
sx q[2];
rz(-1.4497309) q[2];
sx q[2];
rz(2.6775286) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.590893) q[1];
sx q[1];
rz(-1.8807481) q[1];
sx q[1];
rz(2.8551686) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1296785) q[3];
sx q[3];
rz(-2.3647226) q[3];
sx q[3];
rz(-2.7912031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6591349) q[2];
sx q[2];
rz(-1.277801) q[2];
sx q[2];
rz(2.6711312) q[2];
rz(2.2792234) q[3];
sx q[3];
rz(-0.44895288) q[3];
sx q[3];
rz(1.5615777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3622793) q[0];
sx q[0];
rz(-1.1184432) q[0];
sx q[0];
rz(2.733574) q[0];
rz(1.3511924) q[1];
sx q[1];
rz(-1.2974757) q[1];
sx q[1];
rz(-1.8203576) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5274587) q[0];
sx q[0];
rz(-1.2686994) q[0];
sx q[0];
rz(-2.6149349) q[0];
rz(-1.6826019) q[2];
sx q[2];
rz(-2.4833792) q[2];
sx q[2];
rz(-1.3629701) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3776748) q[1];
sx q[1];
rz(-1.498068) q[1];
sx q[1];
rz(-3.0200028) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.42709728) q[3];
sx q[3];
rz(-1.1833041) q[3];
sx q[3];
rz(0.96086335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2230175) q[2];
sx q[2];
rz(-1.7097946) q[2];
sx q[2];
rz(-2.4187386) q[2];
rz(-2.2919848) q[3];
sx q[3];
rz(-1.1725715) q[3];
sx q[3];
rz(-1.0378708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9549114) q[0];
sx q[0];
rz(-2.6246922) q[0];
sx q[0];
rz(0.55996672) q[0];
rz(-2.9365183) q[1];
sx q[1];
rz(-1.2566902) q[1];
sx q[1];
rz(2.8505039) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7322892) q[0];
sx q[0];
rz(-1.8390053) q[0];
sx q[0];
rz(-1.3951673) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.79925691) q[2];
sx q[2];
rz(-1.6148668) q[2];
sx q[2];
rz(1.8792626) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6711802) q[1];
sx q[1];
rz(-0.7268097) q[1];
sx q[1];
rz(-1.8720759) q[1];
rz(-2.970333) q[3];
sx q[3];
rz(-2.3823839) q[3];
sx q[3];
rz(-1.9943135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9186972) q[2];
sx q[2];
rz(-1.6632068) q[2];
sx q[2];
rz(0.29829868) q[2];
rz(-1.1693303) q[3];
sx q[3];
rz(-0.9669286) q[3];
sx q[3];
rz(-1.6303308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(-2.8339612) q[0];
sx q[0];
rz(-1.4606322) q[0];
sx q[0];
rz(2.5391915) q[0];
rz(-2.7700453) q[1];
sx q[1];
rz(-1.5967775) q[1];
sx q[1];
rz(1.0317624) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69799549) q[0];
sx q[0];
rz(-1.5549608) q[0];
sx q[0];
rz(1.3965194) q[0];
rz(-pi) q[1];
rz(-2.6844031) q[2];
sx q[2];
rz(-0.44538272) q[2];
sx q[2];
rz(0.066748652) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6556485) q[1];
sx q[1];
rz(-1.3820501) q[1];
sx q[1];
rz(-0.48635482) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2183275) q[3];
sx q[3];
rz(-2.0482691) q[3];
sx q[3];
rz(1.5058194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3370257) q[2];
sx q[2];
rz(-2.0402543) q[2];
sx q[2];
rz(-1.5186914) q[2];
rz(2.2931781) q[3];
sx q[3];
rz(-2.266326) q[3];
sx q[3];
rz(3.1019822) q[3];
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
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0254211) q[0];
sx q[0];
rz(-0.051055901) q[0];
sx q[0];
rz(2.1436932) q[0];
rz(-2.8669224) q[1];
sx q[1];
rz(-0.88013595) q[1];
sx q[1];
rz(1.7707228) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3034015) q[0];
sx q[0];
rz(-2.5104021) q[0];
sx q[0];
rz(-3.1188008) q[0];
rz(-pi) q[1];
rz(1.9657126) q[2];
sx q[2];
rz(-0.38287258) q[2];
sx q[2];
rz(1.9667786) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.26105598) q[1];
sx q[1];
rz(-0.66268259) q[1];
sx q[1];
rz(0.39775325) q[1];
rz(-0.73565817) q[3];
sx q[3];
rz(-2.1950985) q[3];
sx q[3];
rz(1.5746547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2844598) q[2];
sx q[2];
rz(-2.1993115) q[2];
sx q[2];
rz(3.1119463) q[2];
rz(-2.5263785) q[3];
sx q[3];
rz(-1.7447724) q[3];
sx q[3];
rz(0.342338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089559473) q[0];
sx q[0];
rz(-0.46228662) q[0];
sx q[0];
rz(0.7830559) q[0];
rz(1.1198593) q[1];
sx q[1];
rz(-2.4744787) q[1];
sx q[1];
rz(1.7436183) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8496765) q[0];
sx q[0];
rz(-1.1625942) q[0];
sx q[0];
rz(-2.0565536) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9080284) q[2];
sx q[2];
rz(-1.2937045) q[2];
sx q[2];
rz(1.8943) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.74136342) q[1];
sx q[1];
rz(-1.5796164) q[1];
sx q[1];
rz(2.2386902) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9277738) q[3];
sx q[3];
rz(-0.79584661) q[3];
sx q[3];
rz(-2.4406432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.29584259) q[2];
sx q[2];
rz(-0.58640277) q[2];
sx q[2];
rz(-1.8208549) q[2];
rz(-2.1218421) q[3];
sx q[3];
rz(-1.7176065) q[3];
sx q[3];
rz(1.7216871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19972292) q[0];
sx q[0];
rz(-1.0933192) q[0];
sx q[0];
rz(2.1871908) q[0];
rz(1.154254) q[1];
sx q[1];
rz(-0.64607611) q[1];
sx q[1];
rz(-2.0571713) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53756489) q[0];
sx q[0];
rz(-1.5286235) q[0];
sx q[0];
rz(-2.2794363) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2899288) q[2];
sx q[2];
rz(-1.3251588) q[2];
sx q[2];
rz(-0.40058595) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7560278) q[1];
sx q[1];
rz(-1.4984908) q[1];
sx q[1];
rz(1.0977618) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19799216) q[3];
sx q[3];
rz(-2.3831203) q[3];
sx q[3];
rz(-2.2049979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.083100975) q[2];
sx q[2];
rz(-2.9497171) q[2];
sx q[2];
rz(0.43711883) q[2];
rz(1.7982177) q[3];
sx q[3];
rz(-1.3513887) q[3];
sx q[3];
rz(3.0962211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0992391) q[0];
sx q[0];
rz(-1.0004685) q[0];
sx q[0];
rz(-1.7310671) q[0];
rz(-2.9208185) q[1];
sx q[1];
rz(-0.53123728) q[1];
sx q[1];
rz(-2.208362) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2483053) q[0];
sx q[0];
rz(-0.74549216) q[0];
sx q[0];
rz(3.1214691) q[0];
rz(-pi) q[1];
rz(0.48671203) q[2];
sx q[2];
rz(-0.67811869) q[2];
sx q[2];
rz(-2.202452) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9737967) q[1];
sx q[1];
rz(-0.78404155) q[1];
sx q[1];
rz(1.7586437) q[1];
x q[2];
rz(3.0296928) q[3];
sx q[3];
rz(-0.94116941) q[3];
sx q[3];
rz(-2.1763585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.431388) q[2];
sx q[2];
rz(-0.34222558) q[2];
sx q[2];
rz(0.63146511) q[2];
rz(-2.4225875) q[3];
sx q[3];
rz(-1.7944261) q[3];
sx q[3];
rz(0.24148153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5030293) q[0];
sx q[0];
rz(-1.6521492) q[0];
sx q[0];
rz(-1.8304992) q[0];
rz(0.44166625) q[1];
sx q[1];
rz(-1.8001945) q[1];
sx q[1];
rz(1.1460907) q[1];
rz(1.6497816) q[2];
sx q[2];
rz(-0.38182237) q[2];
sx q[2];
rz(0.67425722) q[2];
rz(-2.5720015) q[3];
sx q[3];
rz(-1.4841595) q[3];
sx q[3];
rz(-1.0356173) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
