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
rz(-0.564841151237488) q[0];
sx q[0];
rz(3.69061270554597) q[0];
sx q[0];
rz(9.24477244018718) q[0];
rz(4.0408673286438) q[1];
sx q[1];
rz(4.00873574812944) q[1];
sx q[1];
rz(8.01981589793369) q[1];
cx q[1],q[0];
rz(4.55668354034424) q[0];
sx q[0];
rz(3.29458326299722) q[0];
sx q[0];
rz(11.5702094793241) q[0];
rz(4.44106674194336) q[2];
sx q[2];
rz(3.71567431290681) q[2];
sx q[2];
rz(11.8070070505063) q[2];
cx q[2],q[1];
rz(-1.63314545154572) q[1];
sx q[1];
rz(4.2594839652353) q[1];
sx q[1];
rz(16.9157323598783) q[1];
rz(3.4698498249054) q[3];
sx q[3];
rz(2.03632882435853) q[3];
sx q[3];
rz(7.80884561537906) q[3];
cx q[3],q[2];
rz(-0.112942419946194) q[2];
sx q[2];
rz(3.88416263659532) q[2];
sx q[2];
rz(13.1479007959287) q[2];
rz(-0.903150081634521) q[3];
sx q[3];
rz(4.79693952401216) q[3];
sx q[3];
rz(9.79188088177844) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.01059317588806) q[0];
sx q[0];
rz(1.40732625325257) q[0];
sx q[0];
rz(11.4291467428128) q[0];
rz(-1.32510232925415) q[1];
sx q[1];
rz(5.61652532418305) q[1];
sx q[1];
rz(10.0890044331472) q[1];
cx q[1],q[0];
rz(0.257608115673065) q[0];
sx q[0];
rz(0.317667873697825) q[0];
sx q[0];
rz(11.5752453565519) q[0];
rz(3.048508644104) q[2];
sx q[2];
rz(0.728903206186839) q[2];
sx q[2];
rz(3.76235959529086) q[2];
cx q[2],q[1];
rz(3.66851186752319) q[1];
sx q[1];
rz(2.01799479325349) q[1];
sx q[1];
rz(12.1651186704557) q[1];
rz(0.31904399394989) q[3];
sx q[3];
rz(3.38226410944993) q[3];
sx q[3];
rz(10.666198706619) q[3];
cx q[3],q[2];
rz(-0.604351580142975) q[2];
sx q[2];
rz(5.12784984906251) q[2];
sx q[2];
rz(10.8058939933698) q[2];
rz(0.857544958591461) q[3];
sx q[3];
rz(0.925679119425364) q[3];
sx q[3];
rz(9.48401363044187) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.183476001024246) q[0];
sx q[0];
rz(7.25586477120454) q[0];
sx q[0];
rz(8.63648281096622) q[0];
rz(3.60860681533813) q[1];
sx q[1];
rz(5.61883440812165) q[1];
sx q[1];
rz(7.12113926409885) q[1];
cx q[1],q[0];
rz(-0.776945114135742) q[0];
sx q[0];
rz(3.5863697548681) q[0];
sx q[0];
rz(10.898866391174) q[0];
rz(0.0639245882630348) q[2];
sx q[2];
rz(5.87104860146577) q[2];
sx q[2];
rz(7.52234003543063) q[2];
cx q[2],q[1];
rz(-4.20875978469849) q[1];
sx q[1];
rz(6.57427898247773) q[1];
sx q[1];
rz(12.8130991220395) q[1];
rz(-1.3784054517746) q[3];
sx q[3];
rz(2.8225401361757) q[3];
sx q[3];
rz(9.72150275706455) q[3];
cx q[3],q[2];
rz(-3.33140110969543) q[2];
sx q[2];
rz(2.80102473695809) q[2];
sx q[2];
rz(11.0298264980237) q[2];
rz(0.324784517288208) q[3];
sx q[3];
rz(1.65805092652375) q[3];
sx q[3];
rz(7.55339334010288) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(3.72454428672791) q[0];
sx q[0];
rz(2.11382094224031) q[0];
sx q[0];
rz(9.47205279245182) q[0];
rz(1.61676967144012) q[1];
sx q[1];
rz(5.84143439133699) q[1];
sx q[1];
rz(7.92225060462161) q[1];
cx q[1],q[0];
rz(1.75280785560608) q[0];
sx q[0];
rz(2.59559169610078) q[0];
sx q[0];
rz(7.30630085467502) q[0];
rz(1.17103981971741) q[2];
sx q[2];
rz(6.77091017563874) q[2];
sx q[2];
rz(9.27356644570037) q[2];
cx q[2],q[1];
rz(-2.66714239120483) q[1];
sx q[1];
rz(6.6339637358957) q[1];
sx q[1];
rz(13.3446120977323) q[1];
rz(-1.03477120399475) q[3];
sx q[3];
rz(3.51754078467424) q[3];
sx q[3];
rz(12.3988291978757) q[3];
cx q[3],q[2];
rz(2.91056847572327) q[2];
sx q[2];
rz(6.67732587655122) q[2];
sx q[2];
rz(13.4768743276517) q[2];
rz(0.856011509895325) q[3];
sx q[3];
rz(1.72497073014314) q[3];
sx q[3];
rz(9.54583967327281) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(4.2471718788147) q[0];
sx q[0];
rz(4.53208306630189) q[0];
sx q[0];
rz(9.35619488953754) q[0];
rz(-0.719064593315125) q[1];
sx q[1];
rz(4.29612782795961) q[1];
sx q[1];
rz(6.70383450984164) q[1];
cx q[1],q[0];
rz(-0.57043582201004) q[0];
sx q[0];
rz(2.12646261056001) q[0];
sx q[0];
rz(8.53464720248386) q[0];
rz(-0.110412076115608) q[2];
sx q[2];
rz(4.93951562245423) q[2];
sx q[2];
rz(10.9874232768933) q[2];
cx q[2],q[1];
rz(-0.334708571434021) q[1];
sx q[1];
rz(5.09839740593965) q[1];
sx q[1];
rz(15.7500447988431) q[1];
rz(1.11334252357483) q[3];
sx q[3];
rz(2.44089994032914) q[3];
sx q[3];
rz(9.82709652780696) q[3];
cx q[3],q[2];
rz(-5.66740226745605) q[2];
sx q[2];
rz(5.31015745003755) q[2];
sx q[2];
rz(8.64512304066821) q[2];
rz(0.0972132161259651) q[3];
sx q[3];
rz(1.32628706296022) q[3];
sx q[3];
rz(10.6212861299436) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.50606107711792) q[0];
sx q[0];
rz(2.47486350138719) q[0];
sx q[0];
rz(8.99789560436412) q[0];
rz(-1.45106685161591) q[1];
sx q[1];
rz(4.26247397263581) q[1];
sx q[1];
rz(6.88765928744479) q[1];
cx q[1],q[0];
rz(2.66478967666626) q[0];
sx q[0];
rz(1.89222970803315) q[0];
sx q[0];
rz(14.1746840238492) q[0];
rz(-3.87046766281128) q[2];
sx q[2];
rz(3.51562512119348) q[2];
sx q[2];
rz(10.6191191434781) q[2];
cx q[2],q[1];
rz(-0.275472432374954) q[1];
sx q[1];
rz(2.94836123486096) q[1];
sx q[1];
rz(6.63493893145725) q[1];
rz(-0.227365538477898) q[3];
sx q[3];
rz(4.16117647488649) q[3];
sx q[3];
rz(8.64944491385623) q[3];
cx q[3],q[2];
rz(0.845885932445526) q[2];
sx q[2];
rz(5.40742579300935) q[2];
sx q[2];
rz(9.17847973703548) q[2];
rz(-1.00364279747009) q[3];
sx q[3];
rz(4.58858147461946) q[3];
sx q[3];
rz(9.34092226474687) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(4.01330709457397) q[0];
sx q[0];
rz(6.35815277894075) q[0];
sx q[0];
rz(12.7713598966519) q[0];
rz(-6.50263261795044) q[1];
sx q[1];
rz(4.58179512818391) q[1];
sx q[1];
rz(12.8457243204038) q[1];
cx q[1],q[0];
rz(1.30592751502991) q[0];
sx q[0];
rz(0.459917457895823) q[0];
sx q[0];
rz(9.68081766962215) q[0];
rz(2.89001965522766) q[2];
sx q[2];
rz(1.61566308339173) q[2];
sx q[2];
rz(10.2363906860273) q[2];
cx q[2],q[1];
rz(3.61259651184082) q[1];
sx q[1];
rz(1.51409724553163) q[1];
sx q[1];
rz(11.5469298124234) q[1];
rz(0.475704818964005) q[3];
sx q[3];
rz(4.98555782635743) q[3];
sx q[3];
rz(9.80862632989093) q[3];
cx q[3],q[2];
rz(2.46519351005554) q[2];
sx q[2];
rz(2.3602843602472) q[2];
sx q[2];
rz(8.54745832680866) q[2];
rz(-1.33897948265076) q[3];
sx q[3];
rz(0.782901438074656) q[3];
sx q[3];
rz(7.67398509978458) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.11246967315674) q[0];
sx q[0];
rz(0.59337297280366) q[0];
sx q[0];
rz(5.94842646121188) q[0];
rz(-2.81074571609497) q[1];
sx q[1];
rz(1.04591623147065) q[1];
sx q[1];
rz(10.0350713491361) q[1];
cx q[1],q[0];
rz(2.42108225822449) q[0];
sx q[0];
rz(1.39255836804444) q[0];
sx q[0];
rz(9.1009454190652) q[0];
rz(-1.91621673107147) q[2];
sx q[2];
rz(4.31717124779756) q[2];
sx q[2];
rz(9.86183775066539) q[2];
cx q[2],q[1];
rz(-5.2122631072998) q[1];
sx q[1];
rz(2.07019320328767) q[1];
sx q[1];
rz(15.4795846700589) q[1];
rz(-1.43707180023193) q[3];
sx q[3];
rz(4.10316559870774) q[3];
sx q[3];
rz(7.80467209815189) q[3];
cx q[3],q[2];
rz(0.941579639911652) q[2];
sx q[2];
rz(4.17850032647187) q[2];
sx q[2];
rz(10.8129915952603) q[2];
rz(0.164182230830193) q[3];
sx q[3];
rz(4.17942682107026) q[3];
sx q[3];
rz(12.2679235696714) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.23090171813965) q[0];
sx q[0];
rz(2.02574840386445) q[0];
sx q[0];
rz(9.28843182920619) q[0];
rz(-3.18960309028625) q[1];
sx q[1];
rz(4.15582826932008) q[1];
sx q[1];
rz(13.8551058530728) q[1];
cx q[1],q[0];
rz(2.89065170288086) q[0];
sx q[0];
rz(4.49552181561524) q[0];
sx q[0];
rz(12.2570130586545) q[0];
rz(-0.451168537139893) q[2];
sx q[2];
rz(4.49942389329011) q[2];
sx q[2];
rz(10.4063962459485) q[2];
cx q[2],q[1];
rz(0.148271888494492) q[1];
sx q[1];
rz(3.43442926009233) q[1];
sx q[1];
rz(10.49222443103) q[1];
rz(2.68194603919983) q[3];
sx q[3];
rz(5.52815213997895) q[3];
sx q[3];
rz(9.98257485627338) q[3];
cx q[3],q[2];
rz(-3.7694149017334) q[2];
sx q[2];
rz(0.00635352929169741) q[2];
sx q[2];
rz(9.60534702836677) q[2];
rz(1.13949656486511) q[3];
sx q[3];
rz(4.65669420559938) q[3];
sx q[3];
rz(12.445990061752) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.9554431438446) q[0];
sx q[0];
rz(5.31016913254792) q[0];
sx q[0];
rz(10.4509395122449) q[0];
rz(1.25632643699646) q[1];
sx q[1];
rz(5.0228740294748) q[1];
sx q[1];
rz(6.21727082728549) q[1];
cx q[1],q[0];
rz(0.724252641201019) q[0];
sx q[0];
rz(5.29588881333406) q[0];
sx q[0];
rz(11.7667062044065) q[0];
rz(-0.171741232275963) q[2];
sx q[2];
rz(4.40318504174287) q[2];
sx q[2];
rz(11.396503186218) q[2];
cx q[2],q[1];
rz(5.89683246612549) q[1];
sx q[1];
rz(0.59659472306306) q[1];
sx q[1];
rz(7.42468140124484) q[1];
rz(0.25392147898674) q[3];
sx q[3];
rz(4.13972774346406) q[3];
sx q[3];
rz(11.9031305074613) q[3];
cx q[3],q[2];
rz(2.57904481887817) q[2];
sx q[2];
rz(2.09360912640626) q[2];
sx q[2];
rz(8.89756486415073) q[2];
rz(0.150500908493996) q[3];
sx q[3];
rz(2.71210572321946) q[3];
sx q[3];
rz(6.63906095027133) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.793138921260834) q[0];
sx q[0];
rz(1.53337994416291) q[0];
sx q[0];
rz(11.7888719797055) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(1.47450315952301) q[1];
sx q[1];
rz(2.27576264937455) q[1];
sx q[1];
rz(11.1977915525357) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-4.52379608154297) q[2];
sx q[2];
rz(6.75806990464265) q[2];
sx q[2];
rz(6.96997973918124) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-0.789535462856293) q[3];
sx q[3];
rz(2.27139315207536) q[3];
sx q[3];
rz(11.1342701673429) q[3];
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
