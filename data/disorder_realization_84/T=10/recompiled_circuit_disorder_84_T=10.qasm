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
rz(0.2660271525383) q[0];
sx q[0];
rz(2.60634830792482) q[0];
sx q[0];
rz(8.6707414150159) q[0];
rz(-5.49298048019409) q[1];
sx q[1];
rz(5.056185515719) q[1];
sx q[1];
rz(8.26391336917087) q[1];
cx q[1],q[0];
rz(2.58295845985413) q[0];
sx q[0];
rz(-1.60038026968902) q[0];
sx q[0];
rz(11.4742810487668) q[0];
rz(3.44731378555298) q[2];
sx q[2];
rz(3.63919753034646) q[2];
sx q[2];
rz(11.3451388835828) q[2];
cx q[2],q[1];
rz(-2.74189925193787) q[1];
sx q[1];
rz(4.20589176018769) q[1];
sx q[1];
rz(7.38110063075229) q[1];
rz(0.824065685272217) q[3];
sx q[3];
rz(6.52517000039155) q[3];
sx q[3];
rz(8.78978208302661) q[3];
cx q[3],q[2];
rz(0.512410879135132) q[2];
sx q[2];
rz(4.43821218808229) q[2];
sx q[2];
rz(10.0912704229276) q[2];
rz(-0.509810566902161) q[3];
sx q[3];
rz(5.13275852997834) q[3];
sx q[3];
rz(10.6982688665311) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.780124127864838) q[0];
sx q[0];
rz(4.54399362404878) q[0];
sx q[0];
rz(8.31208918093845) q[0];
rz(-0.153770461678505) q[1];
sx q[1];
rz(-0.91111168066924) q[1];
sx q[1];
rz(14.3696799039762) q[1];
cx q[1],q[0];
rz(-2.14249444007874) q[0];
sx q[0];
rz(5.55573693116242) q[0];
sx q[0];
rz(7.56323478221103) q[0];
rz(-1.14024639129639) q[2];
sx q[2];
rz(4.02692052920396) q[2];
sx q[2];
rz(11.0986662864606) q[2];
cx q[2],q[1];
rz(6.49191570281982) q[1];
sx q[1];
rz(1.4534229358011) q[1];
sx q[1];
rz(10.6210119485776) q[1];
rz(-1.33093774318695) q[3];
sx q[3];
rz(3.74247506459291) q[3];
sx q[3];
rz(11.6688432454984) q[3];
cx q[3],q[2];
rz(-0.0828369408845901) q[2];
sx q[2];
rz(3.92569431860978) q[2];
sx q[2];
rz(7.46148226260349) q[2];
rz(2.1792049407959) q[3];
sx q[3];
rz(8.3315693457895) q[3];
sx q[3];
rz(13.2318801641385) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.17983794212341) q[0];
sx q[0];
rz(2.63401547272737) q[0];
sx q[0];
rz(7.89774141310855) q[0];
rz(0.642873406410217) q[1];
sx q[1];
rz(4.21811989148194) q[1];
sx q[1];
rz(9.09139270185634) q[1];
cx q[1],q[0];
rz(3.65494179725647) q[0];
sx q[0];
rz(4.91884890397126) q[0];
sx q[0];
rz(10.1242177247922) q[0];
rz(-2.75443911552429) q[2];
sx q[2];
rz(1.69468119938905) q[2];
sx q[2];
rz(7.67141876219913) q[2];
cx q[2],q[1];
rz(0.0458635985851288) q[1];
sx q[1];
rz(3.62155810196931) q[1];
sx q[1];
rz(9.24690658449336) q[1];
rz(0.937366366386414) q[3];
sx q[3];
rz(1.38963499863679) q[3];
sx q[3];
rz(7.34318325518771) q[3];
cx q[3],q[2];
rz(0.12144111096859) q[2];
sx q[2];
rz(4.45554247696931) q[2];
sx q[2];
rz(10.2270946860234) q[2];
rz(-0.241175934672356) q[3];
sx q[3];
rz(5.58435145218904) q[3];
sx q[3];
rz(9.32390703856155) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.08728313446045) q[0];
sx q[0];
rz(4.66489020188386) q[0];
sx q[0];
rz(8.86059353350803) q[0];
rz(3.71972036361694) q[1];
sx q[1];
rz(4.79079047043855) q[1];
sx q[1];
rz(9.93291262387439) q[1];
cx q[1],q[0];
rz(0.807218909263611) q[0];
sx q[0];
rz(2.54265722830827) q[0];
sx q[0];
rz(10.1841440558354) q[0];
rz(-0.932802200317383) q[2];
sx q[2];
rz(4.72984388669068) q[2];
sx q[2];
rz(10.3125390767972) q[2];
cx q[2],q[1];
rz(-0.296246290206909) q[1];
sx q[1];
rz(3.98441913922364) q[1];
sx q[1];
rz(6.63722465037509) q[1];
rz(0.277435064315796) q[3];
sx q[3];
rz(4.53709569771821) q[3];
sx q[3];
rz(9.01207879780933) q[3];
cx q[3],q[2];
rz(1.68867874145508) q[2];
sx q[2];
rz(5.94984355767305) q[2];
sx q[2];
rz(10.0582246541898) q[2];
rz(2.54170346260071) q[3];
sx q[3];
rz(5.13345578511293) q[3];
sx q[3];
rz(10.9250583410184) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.38635814189911) q[0];
sx q[0];
rz(6.58644023736055) q[0];
sx q[0];
rz(9.22867958842918) q[0];
rz(1.26169896125793) q[1];
sx q[1];
rz(3.96207830508287) q[1];
sx q[1];
rz(7.35456321238681) q[1];
cx q[1],q[0];
rz(-2.55398607254028) q[0];
sx q[0];
rz(5.13797214825685) q[0];
sx q[0];
rz(6.96510741709873) q[0];
rz(-0.454829812049866) q[2];
sx q[2];
rz(3.9243393262201) q[2];
sx q[2];
rz(8.96073824762508) q[2];
cx q[2],q[1];
rz(2.56729459762573) q[1];
sx q[1];
rz(4.4291286786371) q[1];
sx q[1];
rz(8.98911879061862) q[1];
rz(5.54049444198608) q[3];
sx q[3];
rz(3.73066672881181) q[3];
sx q[3];
rz(9.03859070538684) q[3];
cx q[3],q[2];
rz(-3.10659575462341) q[2];
sx q[2];
rz(1.80461600621278) q[2];
sx q[2];
rz(10.4067759275357) q[2];
rz(0.185209587216377) q[3];
sx q[3];
rz(5.43920389016206) q[3];
sx q[3];
rz(7.54821727275058) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.53806018829346) q[0];
sx q[0];
rz(4.06101259787614) q[0];
sx q[0];
rz(5.75284264086887) q[0];
rz(-1.84165871143341) q[1];
sx q[1];
rz(4.47136667569215) q[1];
sx q[1];
rz(9.64140520095035) q[1];
cx q[1],q[0];
rz(2.70555710792542) q[0];
sx q[0];
rz(-0.850580541295461) q[0];
sx q[0];
rz(10.6409238338391) q[0];
rz(3.44112706184387) q[2];
sx q[2];
rz(4.10359898407991) q[2];
sx q[2];
rz(8.15760586260959) q[2];
cx q[2],q[1];
rz(-0.129463732242584) q[1];
sx q[1];
rz(6.61858621438081) q[1];
sx q[1];
rz(9.39981919004723) q[1];
rz(-1.37879395484924) q[3];
sx q[3];
rz(-2.2322905937857) q[3];
sx q[3];
rz(13.6830458402555) q[3];
cx q[3],q[2];
rz(2.49143218994141) q[2];
sx q[2];
rz(4.77307191689546) q[2];
sx q[2];
rz(11.5425076246183) q[2];
rz(4.01783084869385) q[3];
sx q[3];
rz(8.72249332268769) q[3];
sx q[3];
rz(8.15321538447543) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.48904705047607) q[0];
sx q[0];
rz(4.30436578591401) q[0];
sx q[0];
rz(8.89584103821918) q[0];
rz(1.61289286613464) q[1];
sx q[1];
rz(4.33382454712922) q[1];
sx q[1];
rz(4.23078439234897) q[1];
cx q[1],q[0];
rz(-2.37503170967102) q[0];
sx q[0];
rz(1.27328172524507) q[0];
sx q[0];
rz(14.143089747421) q[0];
rz(-0.624743342399597) q[2];
sx q[2];
rz(3.95913252432878) q[2];
sx q[2];
rz(11.7848655939023) q[2];
cx q[2],q[1];
rz(0.914280772209167) q[1];
sx q[1];
rz(6.24939742882783) q[1];
sx q[1];
rz(8.04524323939487) q[1];
rz(-1.48005926609039) q[3];
sx q[3];
rz(4.58325055440003) q[3];
sx q[3];
rz(11.1555405616681) q[3];
cx q[3],q[2];
rz(-0.125253245234489) q[2];
sx q[2];
rz(3.94577262003953) q[2];
sx q[2];
rz(7.50870821475192) q[2];
rz(1.6493753194809) q[3];
sx q[3];
rz(4.84632399876649) q[3];
sx q[3];
rz(12.7222451925199) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.656870305538177) q[0];
sx q[0];
rz(6.2211960871988) q[0];
sx q[0];
rz(14.8403382062833) q[0];
rz(3.20881962776184) q[1];
sx q[1];
rz(4.17271641095216) q[1];
sx q[1];
rz(9.22959413229629) q[1];
cx q[1],q[0];
rz(-1.09571003913879) q[0];
sx q[0];
rz(1.36994305451448) q[0];
sx q[0];
rz(11.2958580017011) q[0];
rz(1.60386550426483) q[2];
sx q[2];
rz(5.88184610207612) q[2];
sx q[2];
rz(12.1173121690671) q[2];
cx q[2],q[1];
rz(4.26019668579102) q[1];
sx q[1];
rz(4.67082908947999) q[1];
sx q[1];
rz(8.99491617678806) q[1];
rz(1.92634999752045) q[3];
sx q[3];
rz(2.02046326001222) q[3];
sx q[3];
rz(9.63808763622447) q[3];
cx q[3],q[2];
rz(1.24786698818207) q[2];
sx q[2];
rz(3.35554219980771) q[2];
sx q[2];
rz(13.4718894719998) q[2];
rz(-1.1577821969986) q[3];
sx q[3];
rz(4.31012311776216) q[3];
sx q[3];
rz(9.26562829910918) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.726821959018707) q[0];
sx q[0];
rz(2.04557100136811) q[0];
sx q[0];
rz(11.525418496124) q[0];
rz(3.06293106079102) q[1];
sx q[1];
rz(2.96105957229669) q[1];
sx q[1];
rz(6.63850138186618) q[1];
cx q[1],q[0];
rz(-0.614756345748901) q[0];
sx q[0];
rz(6.66390982468659) q[0];
sx q[0];
rz(9.15846118926212) q[0];
rz(-0.496546566486359) q[2];
sx q[2];
rz(4.92459336121614) q[2];
sx q[2];
rz(10.1324169993322) q[2];
cx q[2],q[1];
rz(4.06027936935425) q[1];
sx q[1];
rz(3.47023800213868) q[1];
sx q[1];
rz(6.1076173543851) q[1];
rz(0.970372140407562) q[3];
sx q[3];
rz(4.55815163453157) q[3];
sx q[3];
rz(7.28359696864291) q[3];
cx q[3],q[2];
rz(0.676552474498749) q[2];
sx q[2];
rz(5.6757380088144) q[2];
sx q[2];
rz(4.87955758570834) q[2];
rz(-3.13535308837891) q[3];
sx q[3];
rz(4.94022861321504) q[3];
sx q[3];
rz(14.1038079023282) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.868260979652405) q[0];
sx q[0];
rz(1.97377589543397) q[0];
sx q[0];
rz(8.12650272845432) q[0];
rz(0.445993304252625) q[1];
sx q[1];
rz(4.25681355794007) q[1];
sx q[1];
rz(6.58392593859836) q[1];
cx q[1],q[0];
rz(2.31076622009277) q[0];
sx q[0];
rz(3.73020801146562) q[0];
sx q[0];
rz(9.00656596421405) q[0];
rz(-0.952063143253326) q[2];
sx q[2];
rz(8.78387227852876) q[2];
sx q[2];
rz(11.9685504197995) q[2];
cx q[2],q[1];
rz(2.44233512878418) q[1];
sx q[1];
rz(4.4211883862787) q[1];
sx q[1];
rz(6.26617000102206) q[1];
rz(2.6273250579834) q[3];
sx q[3];
rz(3.8699798305803) q[3];
sx q[3];
rz(8.15469465254947) q[3];
cx q[3],q[2];
rz(2.51315021514893) q[2];
sx q[2];
rz(2.04378882248933) q[2];
sx q[2];
rz(8.28521630763217) q[2];
rz(1.35097527503967) q[3];
sx q[3];
rz(6.88096013863618) q[3];
sx q[3];
rz(8.25113604067966) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.251912623643875) q[0];
sx q[0];
rz(5.52566924889619) q[0];
sx q[0];
rz(8.92922250031635) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-1.01330137252808) q[1];
sx q[1];
rz(2.69419384201104) q[1];
sx q[1];
rz(13.4631414174955) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(3.19907140731812) q[2];
sx q[2];
rz(5.43711868126924) q[2];
sx q[2];
rz(6.28543255328342) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-0.949284195899963) q[3];
sx q[3];
rz(4.70891455014283) q[3];
sx q[3];
rz(12.3508085966031) q[3];
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
