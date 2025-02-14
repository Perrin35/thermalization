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
rz(1.29901552200317) q[0];
sx q[0];
rz(3.31407797534997) q[0];
sx q[0];
rz(9.40778508073791) q[0];
rz(3.65300822257996) q[1];
sx q[1];
rz(3.63791713316972) q[1];
sx q[1];
rz(12.080892777435) q[1];
cx q[1],q[0];
rz(-0.894285023212433) q[0];
sx q[0];
rz(4.82346955140168) q[0];
sx q[0];
rz(11.1088441371839) q[0];
rz(0.714468598365784) q[2];
sx q[2];
rz(1.03127876122529) q[2];
sx q[2];
rz(4.18175647257968) q[2];
cx q[2],q[1];
rz(-0.790061950683594) q[1];
sx q[1];
rz(2.20215013821656) q[1];
sx q[1];
rz(11.9315576314847) q[1];
rz(1.45620501041412) q[3];
sx q[3];
rz(4.40432360966737) q[3];
sx q[3];
rz(11.3978351116101) q[3];
cx q[3],q[2];
rz(-0.397125065326691) q[2];
sx q[2];
rz(3.03913888533647) q[2];
sx q[2];
rz(8.63219723700687) q[2];
rz(0.986273109912872) q[3];
sx q[3];
rz(7.77055946190888) q[3];
sx q[3];
rz(9.29162473081752) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.886964023113251) q[0];
sx q[0];
rz(1.53810849984223) q[0];
sx q[0];
rz(9.31904055773422) q[0];
rz(3.89950704574585) q[1];
sx q[1];
rz(5.57038298447663) q[1];
sx q[1];
rz(9.90069878696605) q[1];
cx q[1],q[0];
rz(2.02378129959106) q[0];
sx q[0];
rz(6.02896061738069) q[0];
sx q[0];
rz(7.45798573493167) q[0];
rz(-2.9980788230896) q[2];
sx q[2];
rz(5.20102349122102) q[2];
sx q[2];
rz(7.53867170809909) q[2];
cx q[2],q[1];
rz(-1.19131743907928) q[1];
sx q[1];
rz(1.44982674916322) q[1];
sx q[1];
rz(13.1366135835569) q[1];
rz(-0.643893539905548) q[3];
sx q[3];
rz(-1.33126482169097) q[3];
sx q[3];
rz(10.288273549072) q[3];
cx q[3],q[2];
rz(0.146257489919662) q[2];
sx q[2];
rz(6.75759282906587) q[2];
sx q[2];
rz(10.7417080163877) q[2];
rz(2.31381916999817) q[3];
sx q[3];
rz(1.09319487412507) q[3];
sx q[3];
rz(10.2620246171872) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.718831717967987) q[0];
sx q[0];
rz(4.48549285729463) q[0];
sx q[0];
rz(8.18788621424838) q[0];
rz(1.39245665073395) q[1];
sx q[1];
rz(4.862457664805) q[1];
sx q[1];
rz(8.7344444155614) q[1];
cx q[1],q[0];
rz(0.250940978527069) q[0];
sx q[0];
rz(4.25779715378816) q[0];
sx q[0];
rz(10.3713868021886) q[0];
rz(-3.6049108505249) q[2];
sx q[2];
rz(4.70881143410737) q[2];
sx q[2];
rz(7.91332278250858) q[2];
cx q[2],q[1];
rz(2.95430898666382) q[1];
sx q[1];
rz(3.39748230774934) q[1];
sx q[1];
rz(8.56870267390414) q[1];
rz(2.23472571372986) q[3];
sx q[3];
rz(5.43039694626863) q[3];
sx q[3];
rz(14.2552566289823) q[3];
cx q[3],q[2];
rz(3.66007280349731) q[2];
sx q[2];
rz(4.66395023663575) q[2];
sx q[2];
rz(6.20931885241672) q[2];
rz(1.98339235782623) q[3];
sx q[3];
rz(5.06621387799317) q[3];
sx q[3];
rz(13.9733381032865) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.91279196739197) q[0];
sx q[0];
rz(3.08562902559573) q[0];
sx q[0];
rz(6.47611639498874) q[0];
rz(1.8979160785675) q[1];
sx q[1];
rz(4.88697031338746) q[1];
sx q[1];
rz(9.19172864257499) q[1];
cx q[1],q[0];
rz(-3.82843589782715) q[0];
sx q[0];
rz(3.31646134157712) q[0];
sx q[0];
rz(8.42740640639468) q[0];
rz(-0.708791851997375) q[2];
sx q[2];
rz(1.59959116776521) q[2];
sx q[2];
rz(15.1965307950894) q[2];
cx q[2],q[1];
rz(3.45758056640625) q[1];
sx q[1];
rz(4.76866868336732) q[1];
sx q[1];
rz(5.02335164546176) q[1];
rz(4.24165868759155) q[3];
sx q[3];
rz(7.28772226174409) q[3];
sx q[3];
rz(14.4592914342801) q[3];
cx q[3],q[2];
rz(3.35886073112488) q[2];
sx q[2];
rz(4.90381339390809) q[2];
sx q[2];
rz(6.34726903437778) q[2];
rz(3.39281034469604) q[3];
sx q[3];
rz(3.60451129277284) q[3];
sx q[3];
rz(9.71616551875278) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.0837044045329094) q[0];
sx q[0];
rz(4.43523720105226) q[0];
sx q[0];
rz(8.14048669337436) q[0];
rz(3.14902973175049) q[1];
sx q[1];
rz(1.18592718442018) q[1];
sx q[1];
rz(10.381880080692) q[1];
cx q[1],q[0];
rz(1.55146062374115) q[0];
sx q[0];
rz(3.69478681881959) q[0];
sx q[0];
rz(8.78848252295657) q[0];
rz(2.18502330780029) q[2];
sx q[2];
rz(5.22385826905305) q[2];
sx q[2];
rz(7.85865232943698) q[2];
cx q[2],q[1];
rz(-0.196472868323326) q[1];
sx q[1];
rz(5.9396187384897) q[1];
sx q[1];
rz(11.2925965547483) q[1];
rz(1.32940578460693) q[3];
sx q[3];
rz(8.75533881981904) q[3];
sx q[3];
rz(5.4350418805997) q[3];
cx q[3],q[2];
rz(-1.12884223461151) q[2];
sx q[2];
rz(5.57660976250703) q[2];
sx q[2];
rz(7.31350753306552) q[2];
rz(-0.718505918979645) q[3];
sx q[3];
rz(2.05931440194184) q[3];
sx q[3];
rz(10.1313519835393) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.42658233642578) q[0];
sx q[0];
rz(3.88946083386476) q[0];
sx q[0];
rz(9.79188311695262) q[0];
rz(2.78570652008057) q[1];
sx q[1];
rz(5.1658662875467) q[1];
sx q[1];
rz(7.87504122256442) q[1];
cx q[1],q[0];
rz(0.199557602405548) q[0];
sx q[0];
rz(2.39919874270494) q[0];
sx q[0];
rz(10.9385972976606) q[0];
rz(2.9036922454834) q[2];
sx q[2];
rz(3.58937698801095) q[2];
sx q[2];
rz(9.5914583414714) q[2];
cx q[2],q[1];
rz(-1.84252965450287) q[1];
sx q[1];
rz(3.60644391377503) q[1];
sx q[1];
rz(14.6109208822171) q[1];
rz(-1.42450571060181) q[3];
sx q[3];
rz(0.4746357520395) q[3];
sx q[3];
rz(13.2214286088864) q[3];
cx q[3],q[2];
rz(7.25992107391357) q[2];
sx q[2];
rz(4.43473628361756) q[2];
sx q[2];
rz(9.42135037727795) q[2];
rz(2.25547432899475) q[3];
sx q[3];
rz(4.23459997971589) q[3];
sx q[3];
rz(7.72727963923618) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.70573282241821) q[0];
sx q[0];
rz(7.95104184945161) q[0];
sx q[0];
rz(8.70957497357532) q[0];
rz(0.122294783592224) q[1];
sx q[1];
rz(5.2637678702646) q[1];
sx q[1];
rz(9.27715002595588) q[1];
cx q[1],q[0];
rz(1.10811650753021) q[0];
sx q[0];
rz(5.03968611558015) q[0];
sx q[0];
rz(9.78341547249957) q[0];
rz(3.98729205131531) q[2];
sx q[2];
rz(5.72808876832063) q[2];
sx q[2];
rz(8.15702710150882) q[2];
cx q[2],q[1];
rz(1.52772724628448) q[1];
sx q[1];
rz(0.726487787561961) q[1];
sx q[1];
rz(9.70089069604083) q[1];
rz(-0.234695926308632) q[3];
sx q[3];
rz(3.43602615793283) q[3];
sx q[3];
rz(13.460978960983) q[3];
cx q[3],q[2];
rz(2.31323957443237) q[2];
sx q[2];
rz(3.45051896770532) q[2];
sx q[2];
rz(9.53625705688401) q[2];
rz(2.37654232978821) q[3];
sx q[3];
rz(4.56503370602662) q[3];
sx q[3];
rz(12.6002425908963) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.2400803565979) q[0];
sx q[0];
rz(4.11078050931031) q[0];
sx q[0];
rz(11.5635437726895) q[0];
rz(0.785499393939972) q[1];
sx q[1];
rz(4.66648104985292) q[1];
sx q[1];
rz(7.50339601039096) q[1];
cx q[1],q[0];
rz(-1.58831369876862) q[0];
sx q[0];
rz(2.08473411400849) q[0];
sx q[0];
rz(12.5901245832364) q[0];
rz(0.0905916169285774) q[2];
sx q[2];
rz(1.63305095036561) q[2];
sx q[2];
rz(11.1535556077878) q[2];
cx q[2],q[1];
rz(0.739807546138763) q[1];
sx q[1];
rz(0.610817344980784) q[1];
sx q[1];
rz(6.85550401209995) q[1];
rz(2.91369771957397) q[3];
sx q[3];
rz(5.63508764107759) q[3];
sx q[3];
rz(9.57768189012214) q[3];
cx q[3],q[2];
rz(-0.235855236649513) q[2];
sx q[2];
rz(4.88675812085206) q[2];
sx q[2];
rz(9.45857366769716) q[2];
rz(-0.773640513420105) q[3];
sx q[3];
rz(1.61268165906007) q[3];
sx q[3];
rz(13.6291265249173) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.776990234851837) q[0];
sx q[0];
rz(5.66275516350801) q[0];
sx q[0];
rz(9.08239427804157) q[0];
rz(1.42115926742554) q[1];
sx q[1];
rz(5.45992151101167) q[1];
sx q[1];
rz(10.7193209886472) q[1];
cx q[1],q[0];
rz(0.505348205566406) q[0];
sx q[0];
rz(2.74901345570619) q[0];
sx q[0];
rz(7.45920929907962) q[0];
rz(0.849611699581146) q[2];
sx q[2];
rz(7.27685180504853) q[2];
sx q[2];
rz(13.3474704980771) q[2];
cx q[2],q[1];
rz(-2.06599378585815) q[1];
sx q[1];
rz(5.48929754097993) q[1];
sx q[1];
rz(14.3095316648404) q[1];
rz(4.25219440460205) q[3];
sx q[3];
rz(-0.373075334233693) q[3];
sx q[3];
rz(16.4117011785428) q[3];
cx q[3],q[2];
rz(1.29320585727692) q[2];
sx q[2];
rz(4.52698818047578) q[2];
sx q[2];
rz(9.96133652924701) q[2];
rz(3.5544261932373) q[3];
sx q[3];
rz(2.60919133027131) q[3];
sx q[3];
rz(8.31128261088535) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.82630157470703) q[0];
sx q[0];
rz(4.41120520432527) q[0];
sx q[0];
rz(11.069234108917) q[0];
rz(0.810283839702606) q[1];
sx q[1];
rz(4.6981852372461) q[1];
sx q[1];
rz(8.57778981923267) q[1];
cx q[1],q[0];
rz(1.68934559822083) q[0];
sx q[0];
rz(6.04881039460237) q[0];
sx q[0];
rz(8.66518232821628) q[0];
rz(-0.685824275016785) q[2];
sx q[2];
rz(2.36752000649507) q[2];
sx q[2];
rz(14.7333507299344) q[2];
cx q[2],q[1];
rz(0.756397366523743) q[1];
sx q[1];
rz(2.19027307828004) q[1];
sx q[1];
rz(5.95029017924472) q[1];
rz(-1.66693222522736) q[3];
sx q[3];
rz(4.98519006569917) q[3];
sx q[3];
rz(12.0778178930204) q[3];
cx q[3],q[2];
rz(2.10231304168701) q[2];
sx q[2];
rz(1.82590928872163) q[2];
sx q[2];
rz(13.2006041765134) q[2];
rz(2.50411939620972) q[3];
sx q[3];
rz(5.15510741074617) q[3];
sx q[3];
rz(9.20758136212035) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.150443315505981) q[0];
sx q[0];
rz(0.54000964959199) q[0];
sx q[0];
rz(10.2334960460584) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-4.77559757232666) q[1];
sx q[1];
rz(3.50479948719079) q[1];
sx q[1];
rz(11.8881890535276) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(0.45673531293869) q[2];
sx q[2];
rz(4.55788377125794) q[2];
sx q[2];
rz(7.38289711474582) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-0.651675522327423) q[3];
sx q[3];
rz(4.05872878630693) q[3];
sx q[3];
rz(7.5865014552991) q[3];
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
