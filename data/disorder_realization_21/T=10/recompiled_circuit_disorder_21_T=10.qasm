OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.05584597587585) q[0];
sx q[0];
rz(6.36496702035005) q[0];
sx q[0];
rz(9.92624173163577) q[0];
rz(1.49868285655975) q[1];
sx q[1];
rz(3.53775027592713) q[1];
sx q[1];
rz(9.10233936309024) q[1];
cx q[1],q[0];
rz(-1.33813905715942) q[0];
sx q[0];
rz(2.80874663789804) q[0];
sx q[0];
rz(9.68841341733142) q[0];
rz(-2.46039128303528) q[2];
sx q[2];
rz(2.06863644917543) q[2];
sx q[2];
rz(14.3542175054471) q[2];
cx q[2],q[1];
rz(3.44401454925537) q[1];
sx q[1];
rz(4.36612990696961) q[1];
sx q[1];
rz(10.0499025940816) q[1];
rz(0.595789134502411) q[3];
sx q[3];
rz(5.55481305916841) q[3];
sx q[3];
rz(11.1403631925504) q[3];
cx q[3],q[2];
rz(2.63646078109741) q[2];
sx q[2];
rz(2.54870656331117) q[2];
sx q[2];
rz(8.86873821019336) q[2];
rz(2.30891442298889) q[3];
sx q[3];
rz(4.63294640381867) q[3];
sx q[3];
rz(11.6205723047177) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.448220252990723) q[0];
sx q[0];
rz(4.6018592437082) q[0];
sx q[0];
rz(9.58205296694442) q[0];
rz(-2.8804624080658) q[1];
sx q[1];
rz(4.93546059926087) q[1];
sx q[1];
rz(9.53381202965184) q[1];
cx q[1],q[0];
rz(0.635765254497528) q[0];
sx q[0];
rz(6.56942430337007) q[0];
sx q[0];
rz(13.2110969781797) q[0];
rz(0.930958211421967) q[2];
sx q[2];
rz(3.46550905902917) q[2];
sx q[2];
rz(7.77101442813083) q[2];
cx q[2],q[1];
rz(-0.0381395481526852) q[1];
sx q[1];
rz(0.985875757532664) q[1];
sx q[1];
rz(11.9962622880857) q[1];
rz(2.09149575233459) q[3];
sx q[3];
rz(4.56168392499024) q[3];
sx q[3];
rz(9.56129179000064) q[3];
cx q[3],q[2];
rz(1.90335011482239) q[2];
sx q[2];
rz(5.11792567570741) q[2];
sx q[2];
rz(8.16132054328128) q[2];
rz(-0.327109903097153) q[3];
sx q[3];
rz(4.70608309109742) q[3];
sx q[3];
rz(11.3520557641904) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.0644212439656258) q[0];
sx q[0];
rz(3.19088953186805) q[0];
sx q[0];
rz(11.2232393979947) q[0];
rz(-0.247615650296211) q[1];
sx q[1];
rz(5.53654065926606) q[1];
sx q[1];
rz(13.0480441808622) q[1];
cx q[1],q[0];
rz(0.962684810161591) q[0];
sx q[0];
rz(6.91562500794465) q[0];
sx q[0];
rz(11.9045302629392) q[0];
rz(-0.456571251153946) q[2];
sx q[2];
rz(3.80674246151979) q[2];
sx q[2];
rz(11.7262503862302) q[2];
cx q[2],q[1];
rz(-0.183227688074112) q[1];
sx q[1];
rz(5.2540284713083) q[1];
sx q[1];
rz(7.12596032618686) q[1];
rz(0.679694771766663) q[3];
sx q[3];
rz(3.99095234473283) q[3];
sx q[3];
rz(12.0901834726255) q[3];
cx q[3],q[2];
rz(1.80329155921936) q[2];
sx q[2];
rz(7.10058132012422) q[2];
sx q[2];
rz(12.0664787053983) q[2];
rz(2.58061838150024) q[3];
sx q[3];
rz(4.40128806431825) q[3];
sx q[3];
rz(10.885995244972) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.197014763951302) q[0];
sx q[0];
rz(3.30760234792764) q[0];
sx q[0];
rz(8.87259296178027) q[0];
rz(-1.58829700946808) q[1];
sx q[1];
rz(5.38424864609773) q[1];
sx q[1];
rz(7.52792236804172) q[1];
cx q[1],q[0];
rz(-0.23918505012989) q[0];
sx q[0];
rz(3.4810775240236) q[0];
sx q[0];
rz(7.85884044169589) q[0];
rz(-1.44026732444763) q[2];
sx q[2];
rz(1.77032378514344) q[2];
sx q[2];
rz(9.99110493659183) q[2];
cx q[2],q[1];
rz(-0.00754776084795594) q[1];
sx q[1];
rz(4.50050035317475) q[1];
sx q[1];
rz(9.97453562020465) q[1];
rz(1.91361749172211) q[3];
sx q[3];
rz(3.74003401597077) q[3];
sx q[3];
rz(9.53336852639123) q[3];
cx q[3],q[2];
rz(-0.849194705486298) q[2];
sx q[2];
rz(1.25958243210847) q[2];
sx q[2];
rz(11.4156810998838) q[2];
rz(-1.66446471214294) q[3];
sx q[3];
rz(4.7741506417566) q[3];
sx q[3];
rz(9.90772623418971) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.0781590491533279) q[0];
sx q[0];
rz(2.37961909373338) q[0];
sx q[0];
rz(9.50624697505637) q[0];
rz(-3.20405530929565) q[1];
sx q[1];
rz(2.0002580006891) q[1];
sx q[1];
rz(7.92173240183994) q[1];
cx q[1],q[0];
rz(0.0161847285926342) q[0];
sx q[0];
rz(6.88008609612519) q[0];
sx q[0];
rz(12.7949011087339) q[0];
rz(-1.48771107196808) q[2];
sx q[2];
rz(4.82029727299745) q[2];
sx q[2];
rz(8.15807220935031) q[2];
cx q[2],q[1];
rz(1.42457830905914) q[1];
sx q[1];
rz(4.61034479935701) q[1];
sx q[1];
rz(12.8797261476438) q[1];
rz(1.03732657432556) q[3];
sx q[3];
rz(3.71697750886018) q[3];
sx q[3];
rz(7.71716330050632) q[3];
cx q[3],q[2];
rz(0.508267104625702) q[2];
sx q[2];
rz(5.34062901337678) q[2];
sx q[2];
rz(11.3284469604413) q[2];
rz(-1.12269079685211) q[3];
sx q[3];
rz(3.81783065398271) q[3];
sx q[3];
rz(9.94634429215595) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.6102819442749) q[0];
sx q[0];
rz(2.13704815705354) q[0];
sx q[0];
rz(9.69149780868694) q[0];
rz(0.560891270637512) q[1];
sx q[1];
rz(4.98528042634065) q[1];
sx q[1];
rz(11.76783821582) q[1];
cx q[1],q[0];
rz(2.72453665733337) q[0];
sx q[0];
rz(2.8310610969835) q[0];
sx q[0];
rz(9.90101877450153) q[0];
rz(4.21575975418091) q[2];
sx q[2];
rz(4.40271845658357) q[2];
sx q[2];
rz(4.47520253657504) q[2];
cx q[2],q[1];
rz(-3.30679535865784) q[1];
sx q[1];
rz(5.70945731003816) q[1];
sx q[1];
rz(9.52130694537565) q[1];
rz(0.968364238739014) q[3];
sx q[3];
rz(4.03814664681489) q[3];
sx q[3];
rz(9.08112744092151) q[3];
cx q[3],q[2];
rz(-0.160539984703064) q[2];
sx q[2];
rz(5.03421548207337) q[2];
sx q[2];
rz(9.05806139706775) q[2];
rz(-1.26128733158112) q[3];
sx q[3];
rz(4.59492448170716) q[3];
sx q[3];
rz(12.4701630830686) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.409031838178635) q[0];
sx q[0];
rz(4.06186661322648) q[0];
sx q[0];
rz(11.9599835634153) q[0];
rz(0.19730332493782) q[1];
sx q[1];
rz(4.26771810849244) q[1];
sx q[1];
rz(9.88882572054073) q[1];
cx q[1],q[0];
rz(2.79219174385071) q[0];
sx q[0];
rz(2.18786248763139) q[0];
sx q[0];
rz(10.4669952154081) q[0];
rz(-3.22359251976013) q[2];
sx q[2];
rz(3.38317026396329) q[2];
sx q[2];
rz(13.8159708738248) q[2];
cx q[2],q[1];
rz(2.60217666625977) q[1];
sx q[1];
rz(-0.583698598546437) q[1];
sx q[1];
rz(8.85633037089511) q[1];
rz(1.24177348613739) q[3];
sx q[3];
rz(6.03401461442048) q[3];
sx q[3];
rz(13.4065007924955) q[3];
cx q[3],q[2];
rz(-0.934189260005951) q[2];
sx q[2];
rz(5.27999440033967) q[2];
sx q[2];
rz(9.68282333611652) q[2];
rz(1.95596539974213) q[3];
sx q[3];
rz(4.76786819298799) q[3];
sx q[3];
rz(9.42828678562447) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.195140615105629) q[0];
sx q[0];
rz(1.86086812813813) q[0];
sx q[0];
rz(9.04348102807208) q[0];
rz(-0.0952458456158638) q[1];
sx q[1];
rz(4.11402538617189) q[1];
sx q[1];
rz(10.8662957906644) q[1];
cx q[1],q[0];
rz(0.155260741710663) q[0];
sx q[0];
rz(3.07634562452371) q[0];
sx q[0];
rz(8.08772156237766) q[0];
rz(0.958185136318207) q[2];
sx q[2];
rz(4.74379673798616) q[2];
sx q[2];
rz(9.28243134020969) q[2];
cx q[2],q[1];
rz(-3.05635929107666) q[1];
sx q[1];
rz(4.58254483540589) q[1];
sx q[1];
rz(15.8343715429227) q[1];
rz(-1.2368631362915) q[3];
sx q[3];
rz(4.6718294938379) q[3];
sx q[3];
rz(9.55205543934509) q[3];
cx q[3],q[2];
rz(-1.14294457435608) q[2];
sx q[2];
rz(5.87020030816133) q[2];
sx q[2];
rz(9.19818877278968) q[2];
rz(-0.448580265045166) q[3];
sx q[3];
rz(4.74746898015077) q[3];
sx q[3];
rz(10.2545468568723) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.76096928119659) q[0];
sx q[0];
rz(2.38011035521562) q[0];
sx q[0];
rz(8.02577016352817) q[0];
rz(3.45867848396301) q[1];
sx q[1];
rz(4.616683395701) q[1];
sx q[1];
rz(7.26978132723972) q[1];
cx q[1],q[0];
rz(0.875145852565765) q[0];
sx q[0];
rz(3.66189286311204) q[0];
sx q[0];
rz(10.673884010307) q[0];
rz(1.03884887695313) q[2];
sx q[2];
rz(4.08753952582414) q[2];
sx q[2];
rz(6.50553939341708) q[2];
cx q[2],q[1];
rz(-2.03410673141479) q[1];
sx q[1];
rz(1.71205464203889) q[1];
sx q[1];
rz(11.6051268339078) q[1];
rz(3.56241965293884) q[3];
sx q[3];
rz(4.55424872239167) q[3];
sx q[3];
rz(10.9075453042905) q[3];
cx q[3],q[2];
rz(1.92651474475861) q[2];
sx q[2];
rz(3.86859360535676) q[2];
sx q[2];
rz(9.01512131690189) q[2];
rz(0.263272911310196) q[3];
sx q[3];
rz(4.97150150139863) q[3];
sx q[3];
rz(12.0469140767972) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.399193286895752) q[0];
sx q[0];
rz(3.06294585962827) q[0];
sx q[0];
rz(10.8299178838651) q[0];
rz(0.821102261543274) q[1];
sx q[1];
rz(5.36447993119294) q[1];
sx q[1];
rz(8.00927422045871) q[1];
cx q[1],q[0];
rz(3.62564945220947) q[0];
sx q[0];
rz(3.55334744055802) q[0];
sx q[0];
rz(8.47795054911777) q[0];
rz(-3.7299485206604) q[2];
sx q[2];
rz(6.7123214324289) q[2];
sx q[2];
rz(11.6867580175321) q[2];
cx q[2],q[1];
rz(1.07145977020264) q[1];
sx q[1];
rz(1.35766414006288) q[1];
sx q[1];
rz(8.00901589392825) q[1];
rz(2.84709644317627) q[3];
sx q[3];
rz(2.6101971586519) q[3];
sx q[3];
rz(8.94396520256206) q[3];
cx q[3],q[2];
rz(-2.42255091667175) q[2];
sx q[2];
rz(5.98223296006257) q[2];
sx q[2];
rz(9.30067282765313) q[2];
rz(0.965788066387177) q[3];
sx q[3];
rz(4.80870941479737) q[3];
sx q[3];
rz(8.76369521617099) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.32913720607758) q[0];
sx q[0];
rz(1.75826040108735) q[0];
sx q[0];
rz(7.69013748168155) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(3.72091126441956) q[1];
sx q[1];
rz(2.66137993534143) q[1];
sx q[1];
rz(10.1679440498273) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(1.1525970697403) q[2];
sx q[2];
rz(4.54469827015931) q[2];
sx q[2];
rz(8.11884996890231) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(1.60967695713043) q[3];
sx q[3];
rz(1.86741784413392) q[3];
sx q[3];
rz(8.84489921330615) q[3];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
