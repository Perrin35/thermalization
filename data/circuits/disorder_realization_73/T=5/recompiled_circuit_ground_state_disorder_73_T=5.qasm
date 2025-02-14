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
rz(-1.01765275001526) q[0];
sx q[0];
rz(4.05842903454835) q[0];
sx q[0];
rz(9.8596757709901) q[0];
rz(1.59743666648865) q[1];
sx q[1];
rz(3.66060009797151) q[1];
sx q[1];
rz(11.9843477964322) q[1];
cx q[1],q[0];
rz(1.77977418899536) q[0];
sx q[0];
rz(4.86646715004975) q[0];
sx q[0];
rz(9.16071698664829) q[0];
rz(0.596375048160553) q[2];
sx q[2];
rz(2.41084638436372) q[2];
sx q[2];
rz(7.63488707541629) q[2];
cx q[2],q[1];
rz(2.77050447463989) q[1];
sx q[1];
rz(-0.504821149510793) q[1];
sx q[1];
rz(8.89747146367236) q[1];
rz(-1.39277923107147) q[3];
sx q[3];
rz(4.9894689639383) q[3];
sx q[3];
rz(8.15434155463382) q[3];
cx q[3],q[2];
rz(2.58830976486206) q[2];
sx q[2];
rz(4.39901617367799) q[2];
sx q[2];
rz(9.71647671460315) q[2];
rz(-0.450825393199921) q[3];
sx q[3];
rz(2.89016327460343) q[3];
sx q[3];
rz(13.1738161802213) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.701582312583923) q[0];
sx q[0];
rz(5.29271093209321) q[0];
sx q[0];
rz(8.32941994666263) q[0];
rz(-1.19132423400879) q[1];
sx q[1];
rz(4.06159898837144) q[1];
sx q[1];
rz(11.6637947320859) q[1];
cx q[1],q[0];
rz(1.17537367343903) q[0];
sx q[0];
rz(0.928247125940867) q[0];
sx q[0];
rz(10.2187465786855) q[0];
rz(2.77551651000977) q[2];
sx q[2];
rz(0.702090175943919) q[2];
sx q[2];
rz(9.29310607015296) q[2];
cx q[2],q[1];
rz(3.19678592681885) q[1];
sx q[1];
rz(3.45845043857629) q[1];
sx q[1];
rz(10.6003830194394) q[1];
rz(2.07448172569275) q[3];
sx q[3];
rz(4.15183702309663) q[3];
sx q[3];
rz(8.94583991765186) q[3];
cx q[3],q[2];
rz(-0.187874495983124) q[2];
sx q[2];
rz(5.3088847716623) q[2];
sx q[2];
rz(9.54367469846412) q[2];
rz(-0.649057269096375) q[3];
sx q[3];
rz(3.32797414262826) q[3];
sx q[3];
rz(10.8794958352964) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.48979365825653) q[0];
sx q[0];
rz(4.39559498627717) q[0];
sx q[0];
rz(8.90901145934268) q[0];
rz(2.04913973808289) q[1];
sx q[1];
rz(5.87998095353181) q[1];
sx q[1];
rz(11.2911282539289) q[1];
cx q[1],q[0];
rz(2.36088013648987) q[0];
sx q[0];
rz(6.32852593262727) q[0];
sx q[0];
rz(8.26496455668613) q[0];
rz(2.91166806221008) q[2];
sx q[2];
rz(1.87656870682771) q[2];
sx q[2];
rz(8.23354575633212) q[2];
cx q[2],q[1];
rz(-0.66652774810791) q[1];
sx q[1];
rz(4.51629927952821) q[1];
sx q[1];
rz(10.1120635032575) q[1];
rz(0.776712834835052) q[3];
sx q[3];
rz(1.49916723568971) q[3];
sx q[3];
rz(9.01225796937152) q[3];
cx q[3],q[2];
rz(-2.07559967041016) q[2];
sx q[2];
rz(4.45904985268647) q[2];
sx q[2];
rz(11.4065527677457) q[2];
rz(-0.489481508731842) q[3];
sx q[3];
rz(4.69497671921784) q[3];
sx q[3];
rz(10.1445175766866) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.86642611026764) q[0];
sx q[0];
rz(2.81860989530618) q[0];
sx q[0];
rz(9.90643504857227) q[0];
rz(2.21091675758362) q[1];
sx q[1];
rz(4.43845966656739) q[1];
sx q[1];
rz(6.26425144671603) q[1];
cx q[1],q[0];
rz(1.01364946365356) q[0];
sx q[0];
rz(3.19581297610933) q[0];
sx q[0];
rz(9.82176823019191) q[0];
rz(2.47910761833191) q[2];
sx q[2];
rz(4.74728813965852) q[2];
sx q[2];
rz(11.6126363038938) q[2];
cx q[2],q[1];
rz(1.15073096752167) q[1];
sx q[1];
rz(4.12509158452088) q[1];
sx q[1];
rz(10.0562737941663) q[1];
rz(0.732829332351685) q[3];
sx q[3];
rz(4.09380969603593) q[3];
sx q[3];
rz(9.33253587632581) q[3];
cx q[3],q[2];
rz(-0.218617409467697) q[2];
sx q[2];
rz(5.92678180535371) q[2];
sx q[2];
rz(10.1700778365056) q[2];
rz(0.37799146771431) q[3];
sx q[3];
rz(1.99729207356507) q[3];
sx q[3];
rz(10.3133362293164) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.220766305923462) q[0];
sx q[0];
rz(3.4535535295778) q[0];
sx q[0];
rz(10.6156382322232) q[0];
rz(2.65307545661926) q[1];
sx q[1];
rz(2.22564754088456) q[1];
sx q[1];
rz(8.61693576573535) q[1];
cx q[1],q[0];
rz(-0.61208438873291) q[0];
sx q[0];
rz(3.48143911560113) q[0];
sx q[0];
rz(9.14987934230968) q[0];
rz(0.463043361902237) q[2];
sx q[2];
rz(4.50741151173646) q[2];
sx q[2];
rz(9.51740509866878) q[2];
cx q[2],q[1];
rz(2.39920902252197) q[1];
sx q[1];
rz(4.25734809239442) q[1];
sx q[1];
rz(9.32682198881313) q[1];
rz(1.37069892883301) q[3];
sx q[3];
rz(4.09892288048799) q[3];
sx q[3];
rz(8.77028015851184) q[3];
cx q[3],q[2];
rz(-1.00255823135376) q[2];
sx q[2];
rz(5.29046669800813) q[2];
sx q[2];
rz(9.58892883955642) q[2];
rz(0.746033549308777) q[3];
sx q[3];
rz(4.91131428082521) q[3];
sx q[3];
rz(9.49858682452842) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.981107592582703) q[0];
sx q[0];
rz(4.40733400185639) q[0];
sx q[0];
rz(10.1521625280301) q[0];
rz(0.478505939245224) q[1];
sx q[1];
rz(4.59451964695985) q[1];
sx q[1];
rz(10.1616067647855) q[1];
cx q[1],q[0];
rz(1.38723576068878) q[0];
sx q[0];
rz(3.29010172386701) q[0];
sx q[0];
rz(9.80881235598727) q[0];
rz(-0.427992820739746) q[2];
sx q[2];
rz(3.44907212455804) q[2];
sx q[2];
rz(10.2539286971013) q[2];
cx q[2],q[1];
rz(-1.23476541042328) q[1];
sx q[1];
rz(4.3468383868509) q[1];
sx q[1];
rz(11.9117388486783) q[1];
rz(0.387648850679398) q[3];
sx q[3];
rz(4.69949284394319) q[3];
sx q[3];
rz(10.7307371854703) q[3];
cx q[3],q[2];
rz(3.28947067260742) q[2];
sx q[2];
rz(5.25899759133393) q[2];
sx q[2];
rz(10.1099367499272) q[2];
rz(1.00887513160706) q[3];
sx q[3];
rz(5.59312978585298) q[3];
sx q[3];
rz(9.17398501037761) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.20550262928009) q[0];
sx q[0];
rz(6.37714067299897) q[0];
sx q[0];
rz(8.33143994807407) q[0];
rz(3.16534113883972) q[1];
sx q[1];
rz(3.8786268552118) q[1];
sx q[1];
rz(5.78757569789096) q[1];
cx q[1],q[0];
rz(0.281952947378159) q[0];
sx q[0];
rz(6.35034528573091) q[0];
sx q[0];
rz(11.1778193473737) q[0];
rz(1.1497243642807) q[2];
sx q[2];
rz(3.61503252585466) q[2];
sx q[2];
rz(10.4033269643705) q[2];
cx q[2],q[1];
rz(-0.882652819156647) q[1];
sx q[1];
rz(4.53087642987306) q[1];
sx q[1];
rz(10.6409815311353) q[1];
rz(-0.488463133573532) q[3];
sx q[3];
rz(4.32697811921174) q[3];
sx q[3];
rz(12.5559770822446) q[3];
cx q[3],q[2];
rz(0.896398365497589) q[2];
sx q[2];
rz(3.9346828182512) q[2];
sx q[2];
rz(9.12822139858409) q[2];
rz(2.63139295578003) q[3];
sx q[3];
rz(1.52541306813294) q[3];
sx q[3];
rz(9.15335274337932) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.9940013885498) q[0];
sx q[0];
rz(2.18641016085679) q[0];
sx q[0];
rz(10.579091167442) q[0];
rz(1.15559029579163) q[1];
sx q[1];
rz(3.79895201523835) q[1];
sx q[1];
rz(10.8839478254239) q[1];
cx q[1],q[0];
rz(-0.665089130401611) q[0];
sx q[0];
rz(5.5679020007425) q[0];
sx q[0];
rz(11.0305120706479) q[0];
rz(-3.14032125473022) q[2];
sx q[2];
rz(3.86590376694734) q[2];
sx q[2];
rz(11.0437540769498) q[2];
cx q[2],q[1];
rz(-0.369881302118301) q[1];
sx q[1];
rz(5.21272364457185) q[1];
sx q[1];
rz(11.2002106666486) q[1];
rz(0.892919957637787) q[3];
sx q[3];
rz(5.52534118493135) q[3];
sx q[3];
rz(7.37083837985202) q[3];
cx q[3],q[2];
rz(-0.833617210388184) q[2];
sx q[2];
rz(5.74103084405) q[2];
sx q[2];
rz(11.1905271768491) q[2];
rz(0.0863250717520714) q[3];
sx q[3];
rz(3.97426066000993) q[3];
sx q[3];
rz(10.6805699825208) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.19832217693329) q[0];
sx q[0];
rz(3.83540567954118) q[0];
sx q[0];
rz(11.6969182252805) q[0];
rz(-0.729328751564026) q[1];
sx q[1];
rz(3.98516199191148) q[1];
sx q[1];
rz(9.19087643026515) q[1];
cx q[1],q[0];
rz(0.501603543758392) q[0];
sx q[0];
rz(5.86349740822846) q[0];
sx q[0];
rz(11.1129244327466) q[0];
rz(0.0959224700927734) q[2];
sx q[2];
rz(1.36155834992463) q[2];
sx q[2];
rz(12.2659160852353) q[2];
cx q[2],q[1];
rz(1.44713759422302) q[1];
sx q[1];
rz(4.85978797276551) q[1];
sx q[1];
rz(10.5498208761136) q[1];
rz(0.717058300971985) q[3];
sx q[3];
rz(4.95662895043428) q[3];
sx q[3];
rz(7.6007894039075) q[3];
cx q[3],q[2];
rz(-0.0876611545681953) q[2];
sx q[2];
rz(4.78911116917665) q[2];
sx q[2];
rz(8.24853012561008) q[2];
rz(-2.47043991088867) q[3];
sx q[3];
rz(5.03362956841523) q[3];
sx q[3];
rz(7.90867087840244) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.434357076883316) q[0];
sx q[0];
rz(4.003122599917) q[0];
sx q[0];
rz(11.5228295087735) q[0];
rz(3.23889088630676) q[1];
sx q[1];
rz(4.19845178921754) q[1];
sx q[1];
rz(8.30432913302585) q[1];
cx q[1],q[0];
rz(-0.531725108623505) q[0];
sx q[0];
rz(4.77273705800111) q[0];
sx q[0];
rz(10.021051144592) q[0];
rz(2.35666465759277) q[2];
sx q[2];
rz(2.18172696431214) q[2];
sx q[2];
rz(6.4845773935239) q[2];
cx q[2],q[1];
rz(0.0798901170492172) q[1];
sx q[1];
rz(4.88987937768037) q[1];
sx q[1];
rz(9.3752601839523) q[1];
rz(-0.482167780399323) q[3];
sx q[3];
rz(1.40341618855531) q[3];
sx q[3];
rz(9.17901573180362) q[3];
cx q[3],q[2];
rz(1.6843159198761) q[2];
sx q[2];
rz(3.7551094015413) q[2];
sx q[2];
rz(11.2107789277951) q[2];
rz(-1.56543028354645) q[3];
sx q[3];
rz(4.22529080708558) q[3];
sx q[3];
rz(10.4408108949582) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.626513242721558) q[0];
sx q[0];
rz(4.3254201730066) q[0];
sx q[0];
rz(8.26984438895389) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-3.78207540512085) q[1];
sx q[1];
rz(4.20925548871095) q[1];
sx q[1];
rz(8.89747194050952) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(0.522638559341431) q[2];
sx q[2];
rz(3.97319653828675) q[2];
sx q[2];
rz(6.82875726222202) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(0.111079439520836) q[3];
sx q[3];
rz(5.47015801270539) q[3];
sx q[3];
rz(9.16684744357272) q[3];
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
