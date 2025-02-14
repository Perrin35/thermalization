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
rz(0.571036756038666) q[0];
sx q[0];
rz(3.81382259924943) q[0];
sx q[0];
rz(11.0911184310834) q[0];
rz(0.953046798706055) q[1];
sx q[1];
rz(3.3068124969774) q[1];
sx q[1];
rz(10.6934756994168) q[1];
cx q[1],q[0];
rz(0.270862847566605) q[0];
sx q[0];
rz(3.60758796532685) q[0];
sx q[0];
rz(10.4669154643933) q[0];
rz(1.07873272895813) q[2];
sx q[2];
rz(4.11876341898973) q[2];
sx q[2];
rz(10.2037565469663) q[2];
cx q[2],q[1];
rz(-0.424345374107361) q[1];
sx q[1];
rz(4.14589622815187) q[1];
sx q[1];
rz(10.8405924796979) q[1];
rz(1.81469810009003) q[3];
sx q[3];
rz(2.3696405013376) q[3];
sx q[3];
rz(9.41369876488253) q[3];
cx q[3],q[2];
rz(-0.984220087528229) q[2];
sx q[2];
rz(3.48334649403627) q[2];
sx q[2];
rz(10.2179999113004) q[2];
rz(0.673040986061096) q[3];
sx q[3];
rz(4.22701171238954) q[3];
sx q[3];
rz(10.6944551229398) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.157588124275208) q[0];
sx q[0];
rz(3.94420418341691) q[0];
sx q[0];
rz(10.0340870976369) q[0];
rz(-0.209782361984253) q[1];
sx q[1];
rz(3.93292269309098) q[1];
sx q[1];
rz(10.3846638560216) q[1];
cx q[1],q[0];
rz(0.5989089012146) q[0];
sx q[0];
rz(3.23230730195577) q[0];
sx q[0];
rz(10.4821637630384) q[0];
rz(-0.131074965000153) q[2];
sx q[2];
rz(3.63357019622857) q[2];
sx q[2];
rz(10.2080016493718) q[2];
cx q[2],q[1];
rz(-0.349970191717148) q[1];
sx q[1];
rz(3.88117143710191) q[1];
sx q[1];
rz(9.33161748050853) q[1];
rz(1.13502097129822) q[3];
sx q[3];
rz(4.24541989167268) q[3];
sx q[3];
rz(9.34586563556596) q[3];
cx q[3],q[2];
rz(1.35570549964905) q[2];
sx q[2];
rz(3.98087343771989) q[2];
sx q[2];
rz(8.933207756273) q[2];
rz(0.0056886556558311) q[3];
sx q[3];
rz(5.19389692147309) q[3];
sx q[3];
rz(9.93864647149249) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.0950576439499855) q[0];
sx q[0];
rz(3.67704871495301) q[0];
sx q[0];
rz(10.1674127936284) q[0];
rz(0.624780893325806) q[1];
sx q[1];
rz(4.08175251086289) q[1];
sx q[1];
rz(10.4909264802854) q[1];
cx q[1],q[0];
rz(0.103422045707703) q[0];
sx q[0];
rz(3.40684414108331) q[0];
sx q[0];
rz(9.99243203400775) q[0];
rz(-0.320851713418961) q[2];
sx q[2];
rz(4.53143826325471) q[2];
sx q[2];
rz(10.3215422987859) q[2];
cx q[2],q[1];
rz(1.11014926433563) q[1];
sx q[1];
rz(2.71757903893525) q[1];
sx q[1];
rz(9.73480544089481) q[1];
rz(0.393879801034927) q[3];
sx q[3];
rz(2.61632725794847) q[3];
sx q[3];
rz(9.64011200367614) q[3];
cx q[3],q[2];
rz(1.88919579982758) q[2];
sx q[2];
rz(3.33972722490365) q[2];
sx q[2];
rz(8.87220702170535) q[2];
rz(0.507347226142883) q[3];
sx q[3];
rz(4.23078850110108) q[3];
sx q[3];
rz(9.99417886733218) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.660253703594208) q[0];
sx q[0];
rz(3.71891817648942) q[0];
sx q[0];
rz(10.6713551044385) q[0];
rz(1.40568161010742) q[1];
sx q[1];
rz(3.9134853204065) q[1];
sx q[1];
rz(8.70677748917743) q[1];
cx q[1],q[0];
rz(-0.322194904088974) q[0];
sx q[0];
rz(3.9919190128618) q[0];
sx q[0];
rz(10.1246807336728) q[0];
rz(-1.3316605091095) q[2];
sx q[2];
rz(4.16224923928315) q[2];
sx q[2];
rz(10.0118579626004) q[2];
cx q[2],q[1];
rz(0.930546164512634) q[1];
sx q[1];
rz(3.91256103117997) q[1];
sx q[1];
rz(10.7862001418988) q[1];
rz(0.891976833343506) q[3];
sx q[3];
rz(3.34416160185868) q[3];
sx q[3];
rz(9.98442873953983) q[3];
cx q[3],q[2];
rz(1.26728534698486) q[2];
sx q[2];
rz(3.89498070080812) q[2];
sx q[2];
rz(10.12927541732) q[2];
rz(0.382575154304504) q[3];
sx q[3];
rz(4.84512546856935) q[3];
sx q[3];
rz(10.2780271530072) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.480695605278015) q[0];
sx q[0];
rz(3.18179681723053) q[0];
sx q[0];
rz(9.73364139198467) q[0];
rz(2.19246959686279) q[1];
sx q[1];
rz(3.46212026675279) q[1];
sx q[1];
rz(8.87369040249988) q[1];
cx q[1],q[0];
rz(1.23864674568176) q[0];
sx q[0];
rz(4.86469760735566) q[0];
sx q[0];
rz(10.6020660161893) q[0];
rz(1.63961756229401) q[2];
sx q[2];
rz(4.05568739970262) q[2];
sx q[2];
rz(9.06747383474513) q[2];
cx q[2],q[1];
rz(-1.53233623504639) q[1];
sx q[1];
rz(2.92795306642587) q[1];
sx q[1];
rz(11.8103368043821) q[1];
rz(0.794306457042694) q[3];
sx q[3];
rz(4.81563452084596) q[3];
sx q[3];
rz(9.93216792344257) q[3];
cx q[3],q[2];
rz(0.715173125267029) q[2];
sx q[2];
rz(3.89002898533876) q[2];
sx q[2];
rz(10.0435944557111) q[2];
rz(0.623014509677887) q[3];
sx q[3];
rz(3.98351201613481) q[3];
sx q[3];
rz(9.85202965735599) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.0393367856740952) q[0];
sx q[0];
rz(3.80508116086061) q[0];
sx q[0];
rz(9.8866199016492) q[0];
rz(0.496791869401932) q[1];
sx q[1];
rz(4.95962408383424) q[1];
sx q[1];
rz(11.1988225936811) q[1];
cx q[1],q[0];
rz(0.387654662132263) q[0];
sx q[0];
rz(3.82794305880601) q[0];
sx q[0];
rz(9.55026865600749) q[0];
rz(0.276211559772491) q[2];
sx q[2];
rz(2.62819573481614) q[2];
sx q[2];
rz(9.25500289200946) q[2];
cx q[2],q[1];
rz(0.167134717106819) q[1];
sx q[1];
rz(3.44373077352578) q[1];
sx q[1];
rz(9.51883261500999) q[1];
rz(-0.223582610487938) q[3];
sx q[3];
rz(3.77558878262574) q[3];
sx q[3];
rz(9.71828008293315) q[3];
cx q[3],q[2];
rz(-0.10609145462513) q[2];
sx q[2];
rz(4.27365544636781) q[2];
sx q[2];
rz(10.2151481270711) q[2];
rz(0.604493737220764) q[3];
sx q[3];
rz(4.16656711895997) q[3];
sx q[3];
rz(9.84487352370425) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.449643403291702) q[0];
sx q[0];
rz(4.02729430993135) q[0];
sx q[0];
rz(9.8391570508401) q[0];
rz(0.620379328727722) q[1];
sx q[1];
rz(4.36320820649201) q[1];
sx q[1];
rz(7.84199759959384) q[1];
cx q[1],q[0];
rz(0.657114386558533) q[0];
sx q[0];
rz(3.00063909788663) q[0];
sx q[0];
rz(10.3533321380536) q[0];
rz(1.73918235301971) q[2];
sx q[2];
rz(3.03499746521051) q[2];
sx q[2];
rz(8.11838326453372) q[2];
cx q[2],q[1];
rz(1.11902892589569) q[1];
sx q[1];
rz(3.0830585715645) q[1];
sx q[1];
rz(9.26232094167873) q[1];
rz(-0.385635673999786) q[3];
sx q[3];
rz(4.85329976876313) q[3];
sx q[3];
rz(9.69523412584468) q[3];
cx q[3],q[2];
rz(0.15535007417202) q[2];
sx q[2];
rz(3.82939216692979) q[2];
sx q[2];
rz(9.24698466657802) q[2];
rz(0.478846520185471) q[3];
sx q[3];
rz(4.255239995318) q[3];
sx q[3];
rz(9.35647383927509) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.146840453147888) q[0];
sx q[0];
rz(4.12022605736787) q[0];
sx q[0];
rz(9.31220614015266) q[0];
rz(-0.625378906726837) q[1];
sx q[1];
rz(4.26910737355287) q[1];
sx q[1];
rz(10.3140179276387) q[1];
cx q[1],q[0];
rz(0.138769432902336) q[0];
sx q[0];
rz(3.23152415652806) q[0];
sx q[0];
rz(9.92230255006953) q[0];
rz(1.85461008548737) q[2];
sx q[2];
rz(3.41759425600106) q[2];
sx q[2];
rz(10.3641826867978) q[2];
cx q[2],q[1];
rz(1.18046426773071) q[1];
sx q[1];
rz(3.72308442194993) q[1];
sx q[1];
rz(9.86332852243587) q[1];
rz(0.2355887144804) q[3];
sx q[3];
rz(4.77771750290925) q[3];
sx q[3];
rz(10.2470620632093) q[3];
cx q[3],q[2];
rz(1.85829496383667) q[2];
sx q[2];
rz(3.0136643071943) q[2];
sx q[2];
rz(9.69799197315379) q[2];
rz(0.229352712631226) q[3];
sx q[3];
rz(4.72836128075654) q[3];
sx q[3];
rz(10.6933291912) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.447230845689774) q[0];
sx q[0];
rz(2.2706967314058) q[0];
sx q[0];
rz(10.3409930825154) q[0];
rz(-0.182373091578484) q[1];
sx q[1];
rz(3.34781019588048) q[1];
sx q[1];
rz(10.5838165044706) q[1];
cx q[1],q[0];
rz(0.758145093917847) q[0];
sx q[0];
rz(3.44479957421357) q[0];
sx q[0];
rz(10.3111369371335) q[0];
rz(-0.051148809492588) q[2];
sx q[2];
rz(4.61700192292268) q[2];
sx q[2];
rz(9.6408878326337) q[2];
cx q[2],q[1];
rz(0.682848334312439) q[1];
sx q[1];
rz(2.78727513750131) q[1];
sx q[1];
rz(9.89915431141063) q[1];
rz(0.844826221466064) q[3];
sx q[3];
rz(2.07327023346955) q[3];
sx q[3];
rz(9.95564559697315) q[3];
cx q[3],q[2];
rz(0.75659042596817) q[2];
sx q[2];
rz(2.34154102404649) q[2];
sx q[2];
rz(9.10952196120425) q[2];
rz(0.594240188598633) q[3];
sx q[3];
rz(4.02322593529756) q[3];
sx q[3];
rz(10.058880186073) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.8815798163414) q[0];
sx q[0];
rz(3.81027415593202) q[0];
sx q[0];
rz(9.2655385941188) q[0];
rz(1.08819329738617) q[1];
sx q[1];
rz(4.0541075189882) q[1];
sx q[1];
rz(9.69574362634822) q[1];
cx q[1],q[0];
rz(1.2361900806427) q[0];
sx q[0];
rz(4.9632604440027) q[0];
sx q[0];
rz(8.94887433051273) q[0];
rz(1.44225811958313) q[2];
sx q[2];
rz(3.81850543816621) q[2];
sx q[2];
rz(8.51856425999805) q[2];
cx q[2],q[1];
rz(1.38283634185791) q[1];
sx q[1];
rz(3.88462212880189) q[1];
sx q[1];
rz(8.36335680483981) q[1];
rz(0.982175409793854) q[3];
sx q[3];
rz(2.58062198956544) q[3];
sx q[3];
rz(10.0710913300435) q[3];
cx q[3],q[2];
rz(1.34038853645325) q[2];
sx q[2];
rz(3.91375985940034) q[2];
sx q[2];
rz(9.10208050011798) q[2];
rz(0.0580375604331493) q[3];
sx q[3];
rz(3.94312688906724) q[3];
sx q[3];
rz(9.72318536638423) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.829908907413483) q[0];
sx q[0];
rz(3.89822486241395) q[0];
sx q[0];
rz(9.5132926389496) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-0.344874411821365) q[1];
sx q[1];
rz(2.7052145918184) q[1];
sx q[1];
rz(11.0743965864102) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-0.233182430267334) q[2];
sx q[2];
rz(2.86343738635118) q[2];
sx q[2];
rz(10.9450983762662) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(0.0178135931491852) q[3];
sx q[3];
rz(3.79361775715882) q[3];
sx q[3];
rz(11.3376338243405) q[3];
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
