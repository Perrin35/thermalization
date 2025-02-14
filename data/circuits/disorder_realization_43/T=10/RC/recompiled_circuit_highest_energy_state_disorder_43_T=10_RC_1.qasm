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
rz(8.5741841) q[0];
sx q[0];
rz(9.2246715) q[0];
rz(2.9432358) q[1];
sx q[1];
rz(-2.618572) q[1];
sx q[1];
rz(-1.3318292) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1855186) q[0];
sx q[0];
rz(-1.8022493) q[0];
sx q[0];
rz(-1.1197907) q[0];
rz(-pi) q[1];
rz(0.66352377) q[2];
sx q[2];
rz(-0.71394701) q[2];
sx q[2];
rz(2.6235142) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3568452) q[1];
sx q[1];
rz(-0.6797528) q[1];
sx q[1];
rz(-0.076514449) q[1];
rz(-pi) q[2];
x q[2];
rz(0.08272127) q[3];
sx q[3];
rz(-2.045407) q[3];
sx q[3];
rz(3.0623013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.072210463) q[2];
sx q[2];
rz(-0.36082265) q[2];
sx q[2];
rz(1.8270095) q[2];
rz(0.9032816) q[3];
sx q[3];
rz(-1.7777781) q[3];
sx q[3];
rz(-0.34590736) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8697934) q[0];
sx q[0];
rz(-1.485774) q[0];
sx q[0];
rz(-0.85579175) q[0];
rz(-2.3675512) q[1];
sx q[1];
rz(-1.9776521) q[1];
sx q[1];
rz(-2.3537297) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44875328) q[0];
sx q[0];
rz(-1.8215382) q[0];
sx q[0];
rz(-3.1186799) q[0];
rz(-pi) q[1];
rz(-1.6131975) q[2];
sx q[2];
rz(-1.2268492) q[2];
sx q[2];
rz(2.6478772) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1294147) q[1];
sx q[1];
rz(-1.2316634) q[1];
sx q[1];
rz(1.1147092) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66707261) q[3];
sx q[3];
rz(-1.4049585) q[3];
sx q[3];
rz(-1.1950243) q[3];
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
rz(-1.5710477) q[3];
sx q[3];
rz(-1.6377662) q[3];
sx q[3];
rz(1.3554696) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1876672) q[0];
sx q[0];
rz(-1.4952156) q[0];
sx q[0];
rz(-0.060610108) q[0];
rz(-2.4568779) q[1];
sx q[1];
rz(-1.741332) q[1];
sx q[1];
rz(-2.8025467) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6936172) q[0];
sx q[0];
rz(-0.98784271) q[0];
sx q[0];
rz(-1.5956123) q[0];
rz(-1.6940368) q[2];
sx q[2];
rz(-1.3831105) q[2];
sx q[2];
rz(1.0836243) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3587957) q[1];
sx q[1];
rz(-2.7227245) q[1];
sx q[1];
rz(-2.2936506) q[1];
x q[2];
rz(-3.1296785) q[3];
sx q[3];
rz(-0.77687009) q[3];
sx q[3];
rz(2.7912031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6591349) q[2];
sx q[2];
rz(-1.277801) q[2];
sx q[2];
rz(-2.6711312) q[2];
rz(-2.2792234) q[3];
sx q[3];
rz(-0.44895288) q[3];
sx q[3];
rz(1.580015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3622793) q[0];
sx q[0];
rz(-2.0231495) q[0];
sx q[0];
rz(-2.733574) q[0];
rz(-1.7904003) q[1];
sx q[1];
rz(-1.2974757) q[1];
sx q[1];
rz(-1.8203576) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21463284) q[0];
sx q[0];
rz(-2.0713191) q[0];
sx q[0];
rz(1.2248125) q[0];
rz(-1.6826019) q[2];
sx q[2];
rz(-2.4833792) q[2];
sx q[2];
rz(1.7786225) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.76391782) q[1];
sx q[1];
rz(-1.6435247) q[1];
sx q[1];
rz(3.0200028) q[1];
rz(-pi) q[2];
rz(0.42709728) q[3];
sx q[3];
rz(-1.9582886) q[3];
sx q[3];
rz(0.96086335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9185751) q[2];
sx q[2];
rz(-1.431798) q[2];
sx q[2];
rz(-2.4187386) q[2];
rz(0.84960788) q[3];
sx q[3];
rz(-1.9690211) q[3];
sx q[3];
rz(-2.1037219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9549114) q[0];
sx q[0];
rz(-0.51690042) q[0];
sx q[0];
rz(2.5816259) q[0];
rz(2.9365183) q[1];
sx q[1];
rz(-1.2566902) q[1];
sx q[1];
rz(0.29108873) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2084853) q[0];
sx q[0];
rz(-1.7400842) q[0];
sx q[0];
rz(-0.27219682) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.061441378) q[2];
sx q[2];
rz(-2.3413918) q[2];
sx q[2];
rz(-2.7903008) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2771137) q[1];
sx q[1];
rz(-2.2582338) q[1];
sx q[1];
rz(-2.8836125) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75185338) q[3];
sx q[3];
rz(-1.4532148) q[3];
sx q[3];
rz(0.29871179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9186972) q[2];
sx q[2];
rz(-1.6632068) q[2];
sx q[2];
rz(-2.843294) q[2];
rz(-1.9722624) q[3];
sx q[3];
rz(-0.9669286) q[3];
sx q[3];
rz(-1.5112618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8339612) q[0];
sx q[0];
rz(-1.4606322) q[0];
sx q[0];
rz(-2.5391915) q[0];
rz(-2.7700453) q[1];
sx q[1];
rz(-1.5967775) q[1];
sx q[1];
rz(-2.1098302) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1790947) q[0];
sx q[0];
rz(-2.966605) q[0];
sx q[0];
rz(1.4797158) q[0];
x q[1];
rz(1.363109) q[2];
sx q[2];
rz(-1.9676932) q[2];
sx q[2];
rz(0.43235052) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.18373016) q[1];
sx q[1];
rz(-1.0938083) q[1];
sx q[1];
rz(1.7836003) q[1];
x q[2];
rz(-2.5661181) q[3];
sx q[3];
rz(-2.1361975) q[3];
sx q[3];
rz(2.7421212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.80456698) q[2];
sx q[2];
rz(-1.1013384) q[2];
sx q[2];
rz(-1.5186914) q[2];
rz(-0.8484146) q[3];
sx q[3];
rz(-0.87526667) q[3];
sx q[3];
rz(-3.1019822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0254211) q[0];
sx q[0];
rz(-0.051055901) q[0];
sx q[0];
rz(-0.99789944) q[0];
rz(-2.8669224) q[1];
sx q[1];
rz(-2.2614567) q[1];
sx q[1];
rz(-1.7707228) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8557963) q[0];
sx q[0];
rz(-1.5573475) q[0];
sx q[0];
rz(0.6310668) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1758801) q[2];
sx q[2];
rz(-2.7587201) q[2];
sx q[2];
rz(-1.174814) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9129907) q[1];
sx q[1];
rz(-2.1739042) q[1];
sx q[1];
rz(1.8643537) q[1];
rz(-pi) q[2];
rz(-0.79973914) q[3];
sx q[3];
rz(-0.99501401) q[3];
sx q[3];
rz(-2.6588001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.85713282) q[2];
sx q[2];
rz(-2.1993115) q[2];
sx q[2];
rz(-3.1119463) q[2];
rz(2.5263785) q[3];
sx q[3];
rz(-1.7447724) q[3];
sx q[3];
rz(2.7992547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089559473) q[0];
sx q[0];
rz(-0.46228662) q[0];
sx q[0];
rz(-2.3585368) q[0];
rz(1.1198593) q[1];
sx q[1];
rz(-2.4744787) q[1];
sx q[1];
rz(-1.3979744) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8496765) q[0];
sx q[0];
rz(-1.9789985) q[0];
sx q[0];
rz(-1.085039) q[0];
rz(-1.2863761) q[2];
sx q[2];
rz(-1.7952953) q[2];
sx q[2];
rz(-0.38849354) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.82247558) q[1];
sx q[1];
rz(-2.2386595) q[1];
sx q[1];
rz(-0.0112337) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3574202) q[3];
sx q[3];
rz(-2.3436574) q[3];
sx q[3];
rz(-2.7415581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8457501) q[2];
sx q[2];
rz(-0.58640277) q[2];
sx q[2];
rz(-1.3207377) q[2];
rz(-2.1218421) q[3];
sx q[3];
rz(-1.7176065) q[3];
sx q[3];
rz(-1.4199055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9418697) q[0];
sx q[0];
rz(-2.0482735) q[0];
sx q[0];
rz(-0.95440188) q[0];
rz(1.9873387) q[1];
sx q[1];
rz(-0.64607611) q[1];
sx q[1];
rz(2.0571713) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53756489) q[0];
sx q[0];
rz(-1.6129692) q[0];
sx q[0];
rz(0.86215638) q[0];
rz(-2.886359) q[2];
sx q[2];
rz(-1.8430147) q[2];
sx q[2];
rz(-1.24025) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.919405) q[1];
sx q[1];
rz(-2.0424941) q[1];
sx q[1];
rz(-3.0604048) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74862759) q[3];
sx q[3];
rz(-1.4350865) q[3];
sx q[3];
rz(2.3627797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0584917) q[2];
sx q[2];
rz(-2.9497171) q[2];
sx q[2];
rz(-0.43711883) q[2];
rz(-1.3433749) q[3];
sx q[3];
rz(-1.3513887) q[3];
sx q[3];
rz(3.0962211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0992391) q[0];
sx q[0];
rz(-1.0004685) q[0];
sx q[0];
rz(-1.7310671) q[0];
rz(2.9208185) q[1];
sx q[1];
rz(-2.6103554) q[1];
sx q[1];
rz(-2.208362) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8659023) q[0];
sx q[0];
rz(-0.82549107) q[0];
sx q[0];
rz(1.5893713) q[0];
x q[1];
rz(0.61874091) q[2];
sx q[2];
rz(-1.2729984) q[2];
sx q[2];
rz(2.1190475) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.5367914) q[1];
sx q[1];
rz(-1.4385421) q[1];
sx q[1];
rz(-2.3459646) q[1];
rz(1.7228863) q[3];
sx q[3];
rz(-2.503431) q[3];
sx q[3];
rz(-1.1537976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.71020469) q[2];
sx q[2];
rz(-2.7993671) q[2];
sx q[2];
rz(2.5101275) q[2];
rz(2.4225875) q[3];
sx q[3];
rz(-1.3471666) q[3];
sx q[3];
rz(0.24148153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63856335) q[0];
sx q[0];
rz(-1.4894435) q[0];
sx q[0];
rz(1.3110934) q[0];
rz(2.6999264) q[1];
sx q[1];
rz(-1.3413981) q[1];
sx q[1];
rz(-1.9955019) q[1];
rz(-1.491811) q[2];
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
