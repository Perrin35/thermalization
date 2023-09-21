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
rz(0.452861726284027) q[0];
sx q[0];
rz(2.9714669158035) q[0];
sx q[0];
rz(10.2107687354009) q[0];
rz(0.605644702911377) q[1];
sx q[1];
rz(3.79253861506517) q[1];
sx q[1];
rz(8.79069188832446) q[1];
cx q[1],q[0];
rz(0.683610022068024) q[0];
sx q[0];
rz(1.54212203820283) q[0];
sx q[0];
rz(10.1798230171125) q[0];
rz(-1.50793027877808) q[2];
sx q[2];
rz(4.05903688271577) q[2];
sx q[2];
rz(10.6218708515088) q[2];
cx q[2],q[1];
rz(0.699719846248627) q[1];
sx q[1];
rz(2.37082144816453) q[1];
sx q[1];
rz(8.76598641871616) q[1];
rz(1.97098362445831) q[3];
sx q[3];
rz(4.69079068501527) q[3];
sx q[3];
rz(9.24654824136897) q[3];
cx q[3],q[2];
rz(1.22598505020142) q[2];
sx q[2];
rz(5.76610508759553) q[2];
sx q[2];
rz(8.16161391734287) q[2];
rz(1.28492152690887) q[3];
sx q[3];
rz(4.60270980198915) q[3];
sx q[3];
rz(9.13178086876079) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.158535152673721) q[0];
sx q[0];
rz(4.98951950867707) q[0];
sx q[0];
rz(9.86234799622699) q[0];
rz(0.631053984165192) q[1];
sx q[1];
rz(3.47025472124154) q[1];
sx q[1];
rz(12.3452648878019) q[1];
cx q[1],q[0];
rz(0.321999311447144) q[0];
sx q[0];
rz(4.25946608384187) q[0];
sx q[0];
rz(11.2024757623593) q[0];
rz(0.34878870844841) q[2];
sx q[2];
rz(4.81685260136659) q[2];
sx q[2];
rz(10.7054348945539) q[2];
cx q[2],q[1];
rz(-0.0886149555444717) q[1];
sx q[1];
rz(2.66644746263558) q[1];
sx q[1];
rz(11.2631945371549) q[1];
rz(0.097621776163578) q[3];
sx q[3];
rz(2.07115808327729) q[3];
sx q[3];
rz(9.63705476223632) q[3];
cx q[3],q[2];
rz(0.218004032969475) q[2];
sx q[2];
rz(4.57385125954682) q[2];
sx q[2];
rz(10.0130706190984) q[2];
rz(-0.448993921279907) q[3];
sx q[3];
rz(5.85528531868989) q[3];
sx q[3];
rz(11.5183300733487) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.568519711494446) q[0];
sx q[0];
rz(3.04392154713208) q[0];
sx q[0];
rz(10.4252468109052) q[0];
rz(0.725528478622437) q[1];
sx q[1];
rz(4.21613720257814) q[1];
sx q[1];
rz(10.1824768543164) q[1];
cx q[1],q[0];
rz(-1.09627783298492) q[0];
sx q[0];
rz(5.23414650757844) q[0];
sx q[0];
rz(12.1964957475583) q[0];
rz(-1.36331784725189) q[2];
sx q[2];
rz(4.8693198283487) q[2];
sx q[2];
rz(10.0887567758481) q[2];
cx q[2],q[1];
rz(-1.12048995494843) q[1];
sx q[1];
rz(4.24973276455934) q[1];
sx q[1];
rz(11.9089400529782) q[1];
rz(2.20931887626648) q[3];
sx q[3];
rz(4.80200007756288) q[3];
sx q[3];
rz(7.83118221759006) q[3];
cx q[3],q[2];
rz(1.20866012573242) q[2];
sx q[2];
rz(2.91860862274701) q[2];
sx q[2];
rz(11.7299031972806) q[2];
rz(-1.6992734670639) q[3];
sx q[3];
rz(4.61562982399995) q[3];
sx q[3];
rz(9.62860513328716) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.0201241504400969) q[0];
sx q[0];
rz(3.22169802536304) q[0];
sx q[0];
rz(9.93589833974048) q[0];
rz(0.0774436667561531) q[1];
sx q[1];
rz(3.56394607027108) q[1];
sx q[1];
rz(11.2378374099652) q[1];
cx q[1],q[0];
rz(1.11400949954987) q[0];
sx q[0];
rz(2.87675765355165) q[0];
sx q[0];
rz(10.0548221230428) q[0];
rz(-0.494227111339569) q[2];
sx q[2];
rz(5.41246405442292) q[2];
sx q[2];
rz(10.5965182542722) q[2];
cx q[2],q[1];
rz(-1.19917845726013) q[1];
sx q[1];
rz(4.71290257771546) q[1];
sx q[1];
rz(9.15797550081416) q[1];
rz(0.453363686800003) q[3];
sx q[3];
rz(3.38745860953862) q[3];
sx q[3];
rz(9.64072752594157) q[3];
cx q[3],q[2];
rz(-0.0937343090772629) q[2];
sx q[2];
rz(4.35998693306977) q[2];
sx q[2];
rz(8.31771931647464) q[2];
rz(0.472486525774002) q[3];
sx q[3];
rz(4.75622621377046) q[3];
sx q[3];
rz(10.2821528673093) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.0403140187263489) q[0];
sx q[0];
rz(4.0381625016504) q[0];
sx q[0];
rz(9.76297605632945) q[0];
rz(-1.29425537586212) q[1];
sx q[1];
rz(4.7015450318628) q[1];
sx q[1];
rz(9.92928323744937) q[1];
cx q[1],q[0];
rz(0.384066611528397) q[0];
sx q[0];
rz(3.73029354413087) q[0];
sx q[0];
rz(8.99596992730304) q[0];
rz(2.03656601905823) q[2];
sx q[2];
rz(5.05306664307649) q[2];
sx q[2];
rz(9.96094743012592) q[2];
cx q[2],q[1];
rz(-0.0757158473134041) q[1];
sx q[1];
rz(1.26776913006837) q[1];
sx q[1];
rz(10.7456159353177) q[1];
rz(2.23057436943054) q[3];
sx q[3];
rz(4.16447702248628) q[3];
sx q[3];
rz(10.6031090974729) q[3];
cx q[3],q[2];
rz(-2.11110401153564) q[2];
sx q[2];
rz(2.30310735304887) q[2];
sx q[2];
rz(14.0602941274564) q[2];
rz(1.68822872638702) q[3];
sx q[3];
rz(4.07954612572724) q[3];
sx q[3];
rz(10.8070093154828) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.210227385163307) q[0];
sx q[0];
rz(5.40867486794526) q[0];
sx q[0];
rz(10.5156082868497) q[0];
rz(-0.539727210998535) q[1];
sx q[1];
rz(4.29536667664582) q[1];
sx q[1];
rz(9.61356829702064) q[1];
cx q[1],q[0];
rz(0.81297367811203) q[0];
sx q[0];
rz(2.69721737702424) q[0];
sx q[0];
rz(10.3826569080274) q[0];
rz(1.46667921543121) q[2];
sx q[2];
rz(4.44243684609468) q[2];
sx q[2];
rz(7.923452115051) q[2];
cx q[2],q[1];
rz(0.699538826942444) q[1];
sx q[1];
rz(2.890125425654) q[1];
sx q[1];
rz(8.12173352240726) q[1];
rz(2.21419906616211) q[3];
sx q[3];
rz(3.4405484517389) q[3];
sx q[3];
rz(8.26020488738223) q[3];
cx q[3],q[2];
rz(-0.171689003705978) q[2];
sx q[2];
rz(2.11285546620423) q[2];
sx q[2];
rz(7.59831259249851) q[2];
rz(-1.63614892959595) q[3];
sx q[3];
rz(3.49900752504403) q[3];
sx q[3];
rz(11.2830319166104) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.05149066448212) q[0];
sx q[0];
rz(3.97213974793489) q[0];
sx q[0];
rz(9.74877492188617) q[0];
rz(1.84047675132751) q[1];
sx q[1];
rz(5.44787875016267) q[1];
sx q[1];
rz(11.1871777534406) q[1];
cx q[1],q[0];
rz(-2.13490724563599) q[0];
sx q[0];
rz(3.39935964544351) q[0];
sx q[0];
rz(9.75084305404826) q[0];
rz(0.407868087291718) q[2];
sx q[2];
rz(4.2840280850702) q[2];
sx q[2];
rz(10.0637942314069) q[2];
cx q[2],q[1];
rz(0.527232527732849) q[1];
sx q[1];
rz(3.48412466247613) q[1];
sx q[1];
rz(8.35767624377414) q[1];
rz(0.725342452526093) q[3];
sx q[3];
rz(3.92132970889146) q[3];
sx q[3];
rz(9.65985711514159) q[3];
cx q[3],q[2];
rz(0.0327011868357658) q[2];
sx q[2];
rz(4.595935972529) q[2];
sx q[2];
rz(11.6231834649961) q[2];
rz(0.331067830324173) q[3];
sx q[3];
rz(4.91551688511903) q[3];
sx q[3];
rz(9.12965733408138) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.118097640573978) q[0];
sx q[0];
rz(2.55466154416139) q[0];
sx q[0];
rz(10.1890069007795) q[0];
rz(3.28256607055664) q[1];
sx q[1];
rz(3.53766778309877) q[1];
sx q[1];
rz(8.13151512145206) q[1];
cx q[1],q[0];
rz(-0.180235967040062) q[0];
sx q[0];
rz(4.28156009514863) q[0];
sx q[0];
rz(11.213705277435) q[0];
rz(-0.0955351814627647) q[2];
sx q[2];
rz(4.11339953740174) q[2];
sx q[2];
rz(10.1556517839353) q[2];
cx q[2],q[1];
rz(1.58712077140808) q[1];
sx q[1];
rz(1.97381022770936) q[1];
sx q[1];
rz(8.00490734576389) q[1];
rz(-0.200568646192551) q[3];
sx q[3];
rz(4.3200320323282) q[3];
sx q[3];
rz(11.0852882623593) q[3];
cx q[3],q[2];
rz(-0.766888976097107) q[2];
sx q[2];
rz(5.24404135544831) q[2];
sx q[2];
rz(10.0060999751012) q[2];
rz(2.27337217330933) q[3];
sx q[3];
rz(2.72430822451646) q[3];
sx q[3];
rz(7.59461650847598) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.5481036901474) q[0];
sx q[0];
rz(3.72134515841539) q[0];
sx q[0];
rz(10.675421690933) q[0];
rz(-2.14476943016052) q[1];
sx q[1];
rz(4.02528247435624) q[1];
sx q[1];
rz(11.2787579059522) q[1];
cx q[1],q[0];
rz(-1.28519654273987) q[0];
sx q[0];
rz(3.63160604436929) q[0];
sx q[0];
rz(9.47836531921431) q[0];
rz(1.25160813331604) q[2];
sx q[2];
rz(3.64938548405702) q[2];
sx q[2];
rz(9.57427251934215) q[2];
cx q[2],q[1];
rz(0.950836420059204) q[1];
sx q[1];
rz(4.50538554986055) q[1];
sx q[1];
rz(11.914740061752) q[1];
rz(-0.746315062046051) q[3];
sx q[3];
rz(4.67356029351289) q[3];
sx q[3];
rz(8.5908287525098) q[3];
cx q[3],q[2];
rz(-1.32049703598022) q[2];
sx q[2];
rz(3.52950766881044) q[2];
sx q[2];
rz(11.3098956108014) q[2];
rz(0.16658940911293) q[3];
sx q[3];
rz(4.71825972397859) q[3];
sx q[3];
rz(10.4973497152249) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.880841076374054) q[0];
sx q[0];
rz(3.94576284487779) q[0];
sx q[0];
rz(9.74999750255748) q[0];
rz(-1.13517689704895) q[1];
sx q[1];
rz(2.46446410019929) q[1];
sx q[1];
rz(9.05759084819957) q[1];
cx q[1],q[0];
rz(-0.678839147090912) q[0];
sx q[0];
rz(3.19879227330024) q[0];
sx q[0];
rz(10.4202349543492) q[0];
rz(2.19796800613403) q[2];
sx q[2];
rz(4.48331287701661) q[2];
sx q[2];
rz(10.4939961194913) q[2];
cx q[2],q[1];
rz(-2.23284077644348) q[1];
sx q[1];
rz(2.45262608130509) q[1];
sx q[1];
rz(10.6784367322843) q[1];
rz(0.641310393810272) q[3];
sx q[3];
rz(4.31635680993135) q[3];
sx q[3];
rz(9.5158553481023) q[3];
cx q[3],q[2];
rz(0.245115891098976) q[2];
sx q[2];
rz(3.95066175063188) q[2];
sx q[2];
rz(7.34932038783237) q[2];
rz(0.0792447626590729) q[3];
sx q[3];
rz(3.74437311490113) q[3];
sx q[3];
rz(7.94848833083316) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.0358945950865746) q[0];
sx q[0];
rz(3.12512363133068) q[0];
sx q[0];
rz(10.7362776756208) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(2.31300067901611) q[1];
sx q[1];
rz(1.28352204163606) q[1];
sx q[1];
rz(7.31579349040195) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-0.775002717971802) q[2];
sx q[2];
rz(3.19967549865181) q[2];
sx q[2];
rz(7.80769464968845) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-0.704976975917816) q[3];
sx q[3];
rz(2.22720745404298) q[3];
sx q[3];
rz(11.3401039600293) q[3];
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