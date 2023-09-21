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
rz(-0.335586935281754) q[0];
sx q[0];
rz(4.08685627778108) q[0];
sx q[0];
rz(9.95036903618976) q[0];
rz(0.243191495537758) q[1];
sx q[1];
rz(4.37421241601045) q[1];
sx q[1];
rz(10.3296196818273) q[1];
cx q[1],q[0];
rz(-0.685695052146912) q[0];
sx q[0];
rz(4.21822431881959) q[0];
sx q[0];
rz(10.5003454446714) q[0];
rz(2.28440427780151) q[2];
sx q[2];
rz(3.43701595266397) q[2];
sx q[2];
rz(7.67403886317416) q[2];
cx q[2],q[1];
rz(0.799965858459473) q[1];
sx q[1];
rz(3.60871881444985) q[1];
sx q[1];
rz(8.881005203716) q[1];
rz(-0.537006556987762) q[3];
sx q[3];
rz(3.4992423077398) q[3];
sx q[3];
rz(10.7305528879087) q[3];
cx q[3],q[2];
rz(0.20377542078495) q[2];
sx q[2];
rz(4.87693873246247) q[2];
sx q[2];
rz(9.52102233319684) q[2];
rz(1.03592669963837) q[3];
sx q[3];
rz(3.52873552043969) q[3];
sx q[3];
rz(9.57849214076206) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.440898329019547) q[0];
sx q[0];
rz(3.53273290594155) q[0];
sx q[0];
rz(10.1899491310041) q[0];
rz(-1.84934294223785) q[1];
sx q[1];
rz(3.62680244644219) q[1];
sx q[1];
rz(11.9034204244535) q[1];
cx q[1],q[0];
rz(-0.637558579444885) q[0];
sx q[0];
rz(3.23248163064057) q[0];
sx q[0];
rz(10.1779333114545) q[0];
rz(1.16947150230408) q[2];
sx q[2];
rz(5.3255571444803) q[2];
sx q[2];
rz(10.2319078206937) q[2];
cx q[2],q[1];
rz(-0.0291649736464024) q[1];
sx q[1];
rz(4.79487613041932) q[1];
sx q[1];
rz(8.2957561969678) q[1];
rz(0.161969602108002) q[3];
sx q[3];
rz(4.23099079926545) q[3];
sx q[3];
rz(9.37932591362997) q[3];
cx q[3],q[2];
rz(-0.0489166527986526) q[2];
sx q[2];
rz(4.19016304810578) q[2];
sx q[2];
rz(9.09576398729488) q[2];
rz(-0.665502309799194) q[3];
sx q[3];
rz(2.92329590220983) q[3];
sx q[3];
rz(10.7425198316495) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.592778503894806) q[0];
sx q[0];
rz(2.91872616310651) q[0];
sx q[0];
rz(9.64711812733814) q[0];
rz(1.01733434200287) q[1];
sx q[1];
rz(3.86287388403947) q[1];
sx q[1];
rz(9.94346580504581) q[1];
cx q[1],q[0];
rz(1.40435361862183) q[0];
sx q[0];
rz(3.95599207480485) q[0];
sx q[0];
rz(9.92435354589626) q[0];
rz(1.90962946414948) q[2];
sx q[2];
rz(3.42300126154954) q[2];
sx q[2];
rz(7.40103743075534) q[2];
cx q[2],q[1];
rz(2.05734038352966) q[1];
sx q[1];
rz(3.94735393126542) q[1];
sx q[1];
rz(9.18436401187583) q[1];
rz(1.58549892902374) q[3];
sx q[3];
rz(3.19621010695631) q[3];
sx q[3];
rz(10.285225725166) q[3];
cx q[3],q[2];
rz(1.28059327602386) q[2];
sx q[2];
rz(3.50342062314088) q[2];
sx q[2];
rz(9.35623442976877) q[2];
rz(0.602442562580109) q[3];
sx q[3];
rz(2.37903627951676) q[3];
sx q[3];
rz(9.56379712223216) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.49087262153625) q[0];
sx q[0];
rz(3.91871443589265) q[0];
sx q[0];
rz(9.25053288637801) q[0];
rz(0.530255913734436) q[1];
sx q[1];
rz(4.62418022950227) q[1];
sx q[1];
rz(9.93786898850604) q[1];
cx q[1],q[0];
rz(0.952723205089569) q[0];
sx q[0];
rz(5.33085480530793) q[0];
sx q[0];
rz(9.52185314743921) q[0];
rz(0.763637959957123) q[2];
sx q[2];
rz(4.24611154397065) q[2];
sx q[2];
rz(10.9640567064206) q[2];
cx q[2],q[1];
rz(0.404011338949203) q[1];
sx q[1];
rz(2.02391317685182) q[1];
sx q[1];
rz(9.27092360555335) q[1];
rz(0.206391468644142) q[3];
sx q[3];
rz(4.74257400830323) q[3];
sx q[3];
rz(9.57489741443797) q[3];
cx q[3],q[2];
rz(2.66673684120178) q[2];
sx q[2];
rz(3.50777092774446) q[2];
sx q[2];
rz(9.65466493963405) q[2];
rz(-0.419044673442841) q[3];
sx q[3];
rz(4.4906326850229) q[3];
sx q[3];
rz(9.88405594824954) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.672236502170563) q[0];
sx q[0];
rz(3.89278897841508) q[0];
sx q[0];
rz(10.1015972256581) q[0];
rz(0.493044018745422) q[1];
sx q[1];
rz(4.09058425028855) q[1];
sx q[1];
rz(10.040838277332) q[1];
cx q[1],q[0];
rz(0.842329323291779) q[0];
sx q[0];
rz(3.20512517740066) q[0];
sx q[0];
rz(10.024175977699) q[0];
rz(1.26625394821167) q[2];
sx q[2];
rz(3.52282309730584) q[2];
sx q[2];
rz(9.19064096211597) q[2];
cx q[2],q[1];
rz(0.383560091257095) q[1];
sx q[1];
rz(3.7192349751764) q[1];
sx q[1];
rz(11.0056704044263) q[1];
rz(0.676643550395966) q[3];
sx q[3];
rz(4.47072330315644) q[3];
sx q[3];
rz(11.0040709733884) q[3];
cx q[3],q[2];
rz(0.734135091304779) q[2];
sx q[2];
rz(2.58717581828172) q[2];
sx q[2];
rz(9.16941068171664) q[2];
rz(1.53639614582062) q[3];
sx q[3];
rz(4.09401038487489) q[3];
sx q[3];
rz(8.65068171023532) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.422463208436966) q[0];
sx q[0];
rz(4.10779038270051) q[0];
sx q[0];
rz(9.897278225414) q[0];
rz(0.526081442832947) q[1];
sx q[1];
rz(3.35145062406594) q[1];
sx q[1];
rz(11.6816038846891) q[1];
cx q[1],q[0];
rz(-0.173818290233612) q[0];
sx q[0];
rz(3.52982980211312) q[0];
sx q[0];
rz(9.62149957417651) q[0];
rz(0.014147368259728) q[2];
sx q[2];
rz(4.50235119660432) q[2];
sx q[2];
rz(10.5673916101377) q[2];
cx q[2],q[1];
rz(0.519845187664032) q[1];
sx q[1];
rz(4.42487707932527) q[1];
sx q[1];
rz(9.85255709885761) q[1];
rz(-0.401745527982712) q[3];
sx q[3];
rz(2.57079765399034) q[3];
sx q[3];
rz(10.8874162197034) q[3];
cx q[3],q[2];
rz(-0.744495987892151) q[2];
sx q[2];
rz(3.94656625588471) q[2];
sx q[2];
rz(9.73615626095935) q[2];
rz(1.77296757698059) q[3];
sx q[3];
rz(2.68406483729417) q[3];
sx q[3];
rz(9.95409826039478) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.418065279722214) q[0];
sx q[0];
rz(2.86302456458146) q[0];
sx q[0];
rz(9.48584861158534) q[0];
rz(0.0401849895715714) q[1];
sx q[1];
rz(4.30270007451112) q[1];
sx q[1];
rz(10.1576704144399) q[1];
cx q[1],q[0];
rz(0.311567634344101) q[0];
sx q[0];
rz(3.82810738881166) q[0];
sx q[0];
rz(10.3410699725072) q[0];
rz(-0.549675524234772) q[2];
sx q[2];
rz(3.92293086846406) q[2];
sx q[2];
rz(9.90652812122508) q[2];
cx q[2],q[1];
rz(1.14547538757324) q[1];
sx q[1];
rz(3.76737591822679) q[1];
sx q[1];
rz(9.41362339760318) q[1];
rz(0.366883665323257) q[3];
sx q[3];
rz(3.25652737368877) q[3];
sx q[3];
rz(10.1113298296849) q[3];
cx q[3],q[2];
rz(1.07330334186554) q[2];
sx q[2];
rz(4.33855191071565) q[2];
sx q[2];
rz(9.93936160802051) q[2];
rz(1.23752808570862) q[3];
sx q[3];
rz(3.72461095650727) q[3];
sx q[3];
rz(8.87985960244342) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.0562036894261837) q[0];
sx q[0];
rz(4.79789284070069) q[0];
sx q[0];
rz(10.1342046618383) q[0];
rz(1.50526940822601) q[1];
sx q[1];
rz(4.21535232861573) q[1];
sx q[1];
rz(9.7034910082738) q[1];
cx q[1],q[0];
rz(1.3942791223526) q[0];
sx q[0];
rz(2.89076674182946) q[0];
sx q[0];
rz(10.4210626840512) q[0];
rz(1.01217222213745) q[2];
sx q[2];
rz(4.92015627224977) q[2];
sx q[2];
rz(9.07697588800594) q[2];
cx q[2],q[1];
rz(1.7396000623703) q[1];
sx q[1];
rz(3.52471637924249) q[1];
sx q[1];
rz(8.8086039185445) q[1];
rz(0.33155357837677) q[3];
sx q[3];
rz(3.37729773123796) q[3];
sx q[3];
rz(10.6577447414319) q[3];
cx q[3],q[2];
rz(-0.301482081413269) q[2];
sx q[2];
rz(4.43953362305696) q[2];
sx q[2];
rz(9.48085745646759) q[2];
rz(0.855148315429688) q[3];
sx q[3];
rz(3.59007546504075) q[3];
sx q[3];
rz(9.82996747492954) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.994960188865662) q[0];
sx q[0];
rz(3.31790060003335) q[0];
sx q[0];
rz(10.3936707734983) q[0];
rz(-0.473371982574463) q[1];
sx q[1];
rz(3.9362087567621) q[1];
sx q[1];
rz(10.5673125743787) q[1];
cx q[1],q[0];
rz(0.245070949196815) q[0];
sx q[0];
rz(3.45043570001657) q[0];
sx q[0];
rz(9.72772443889781) q[0];
rz(1.57368659973145) q[2];
sx q[2];
rz(3.84493872721727) q[2];
sx q[2];
rz(9.32187685220643) q[2];
cx q[2],q[1];
rz(-1.68563187122345) q[1];
sx q[1];
rz(3.69614115555818) q[1];
sx q[1];
rz(11.2355149745862) q[1];
rz(0.254083245992661) q[3];
sx q[3];
rz(3.82772466738755) q[3];
sx q[3];
rz(9.80954394339725) q[3];
cx q[3],q[2];
rz(1.24509763717651) q[2];
sx q[2];
rz(1.92180398304994) q[2];
sx q[2];
rz(10.4455600738446) q[2];
rz(0.323723703622818) q[3];
sx q[3];
rz(3.89457836945588) q[3];
sx q[3];
rz(9.43595007564082) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.979941189289093) q[0];
sx q[0];
rz(3.16949088883633) q[0];
sx q[0];
rz(10.1261836647908) q[0];
rz(0.91570633649826) q[1];
sx q[1];
rz(5.2748750766092) q[1];
sx q[1];
rz(8.18621537684604) q[1];
cx q[1],q[0];
rz(0.21047131717205) q[0];
sx q[0];
rz(3.70874807436998) q[0];
sx q[0];
rz(9.83428833483859) q[0];
rz(-0.0513259992003441) q[2];
sx q[2];
rz(3.96005544264848) q[2];
sx q[2];
rz(10.2135749816815) q[2];
cx q[2],q[1];
rz(0.306759476661682) q[1];
sx q[1];
rz(3.80814162095124) q[1];
sx q[1];
rz(9.23381324707671) q[1];
rz(-0.168335303664207) q[3];
sx q[3];
rz(3.7406891306215) q[3];
sx q[3];
rz(10.825569009773) q[3];
cx q[3],q[2];
rz(0.686765849590302) q[2];
sx q[2];
rz(4.23802378972108) q[2];
sx q[2];
rz(9.38093942253991) q[2];
rz(1.20101261138916) q[3];
sx q[3];
rz(3.87692973216111) q[3];
sx q[3];
rz(10.428362941734) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.94168496131897) q[0];
sx q[0];
rz(4.00481203396852) q[0];
sx q[0];
rz(10.816607451431) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(0.844746947288513) q[1];
sx q[1];
rz(2.74389103253419) q[1];
sx q[1];
rz(10.1657372474591) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(0.961610853672028) q[2];
sx q[2];
rz(4.08939084609086) q[2];
sx q[2];
rz(9.77514249681636) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(1.23522102832794) q[3];
sx q[3];
rz(2.85762909253175) q[3];
sx q[3];
rz(8.97275022267505) q[3];
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
