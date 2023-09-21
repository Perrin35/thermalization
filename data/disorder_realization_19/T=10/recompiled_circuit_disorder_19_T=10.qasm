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
rz(1.02695453166962) q[0];
sx q[0];
rz(4.59335246880586) q[0];
sx q[0];
rz(10.0703548550527) q[0];
rz(0.378807008266449) q[1];
sx q[1];
rz(4.91034761269624) q[1];
sx q[1];
rz(11.0684392213742) q[1];
cx q[1],q[0];
rz(-0.492555320262909) q[0];
sx q[0];
rz(4.85678509076173) q[0];
sx q[0];
rz(9.06494966744586) q[0];
rz(-0.0954043567180634) q[2];
sx q[2];
rz(5.65049782593782) q[2];
sx q[2];
rz(9.64581890999481) q[2];
cx q[2],q[1];
rz(-2.54219961166382) q[1];
sx q[1];
rz(4.92717984517152) q[1];
sx q[1];
rz(13.3944866418759) q[1];
rz(0.399980485439301) q[3];
sx q[3];
rz(4.10225638945634) q[3];
sx q[3];
rz(8.74893251656696) q[3];
cx q[3],q[2];
rz(2.22151851654053) q[2];
sx q[2];
rz(4.66207710106904) q[2];
sx q[2];
rz(12.9067795038144) q[2];
rz(-0.832996249198914) q[3];
sx q[3];
rz(0.703279884653636) q[3];
sx q[3];
rz(11.2097060441892) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.672656357288361) q[0];
sx q[0];
rz(2.75149727066094) q[0];
sx q[0];
rz(9.96976385115787) q[0];
rz(0.908220589160919) q[1];
sx q[1];
rz(3.46908632119233) q[1];
sx q[1];
rz(10.2497370004575) q[1];
cx q[1],q[0];
rz(-1.41754066944122) q[0];
sx q[0];
rz(5.11528971989686) q[0];
sx q[0];
rz(7.65035519599124) q[0];
rz(1.97370147705078) q[2];
sx q[2];
rz(4.7219719012552) q[2];
sx q[2];
rz(9.11068896054431) q[2];
cx q[2],q[1];
rz(2.0032012462616) q[1];
sx q[1];
rz(1.01112857659394) q[1];
sx q[1];
rz(9.5619155973117) q[1];
rz(1.7663232088089) q[3];
sx q[3];
rz(4.15447738965089) q[3];
sx q[3];
rz(10.0248025417249) q[3];
cx q[3],q[2];
rz(2.55435824394226) q[2];
sx q[2];
rz(3.6645583828264) q[2];
sx q[2];
rz(6.10677287577792) q[2];
rz(-0.130886256694794) q[3];
sx q[3];
rz(4.31708052952821) q[3];
sx q[3];
rz(9.45743524505898) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.7628870010376) q[0];
sx q[0];
rz(2.55809065897996) q[0];
sx q[0];
rz(8.95503080486461) q[0];
rz(-1.52479553222656) q[1];
sx q[1];
rz(2.66416147549684) q[1];
sx q[1];
rz(12.5278310537259) q[1];
cx q[1],q[0];
rz(2.8124213218689) q[0];
sx q[0];
rz(3.1338131084391) q[0];
sx q[0];
rz(9.50727574377462) q[0];
rz(-0.153624624013901) q[2];
sx q[2];
rz(4.9405678828531) q[2];
sx q[2];
rz(8.01197538375064) q[2];
cx q[2],q[1];
rz(0.0270278342068195) q[1];
sx q[1];
rz(4.77389696438844) q[1];
sx q[1];
rz(9.91485652922794) q[1];
rz(0.774499416351318) q[3];
sx q[3];
rz(2.37818035681779) q[3];
sx q[3];
rz(11.6754939317624) q[3];
cx q[3],q[2];
rz(3.72435021400452) q[2];
sx q[2];
rz(4.23113385041291) q[2];
sx q[2];
rz(10.3373020052831) q[2];
rz(1.83306658267975) q[3];
sx q[3];
rz(5.2794670184427) q[3];
sx q[3];
rz(10.8661668062131) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.383476048707962) q[0];
sx q[0];
rz(4.45706537564332) q[0];
sx q[0];
rz(9.77447569965526) q[0];
rz(-1.89674592018127) q[1];
sx q[1];
rz(3.44520774682099) q[1];
sx q[1];
rz(12.7615294217984) q[1];
cx q[1],q[0];
rz(0.124766066670418) q[0];
sx q[0];
rz(3.29365730484063) q[0];
sx q[0];
rz(11.9102329969327) q[0];
rz(2.05173134803772) q[2];
sx q[2];
rz(4.22009685833985) q[2];
sx q[2];
rz(8.99875105022594) q[2];
cx q[2],q[1];
rz(0.819820404052734) q[1];
sx q[1];
rz(1.12045493920381) q[1];
sx q[1];
rz(8.11576697825595) q[1];
rz(0.016276553273201) q[3];
sx q[3];
rz(4.07332590420777) q[3];
sx q[3];
rz(8.40134522914096) q[3];
cx q[3],q[2];
rz(-0.230947092175484) q[2];
sx q[2];
rz(5.79329100449617) q[2];
sx q[2];
rz(13.8341135740201) q[2];
rz(-1.06864821910858) q[3];
sx q[3];
rz(4.15410772164399) q[3];
sx q[3];
rz(10.2651141643445) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.1068263053894) q[0];
sx q[0];
rz(4.59383955796296) q[0];
sx q[0];
rz(12.4267573118131) q[0];
rz(-1.07682490348816) q[1];
sx q[1];
rz(2.10075167019899) q[1];
sx q[1];
rz(12.1882793664853) q[1];
cx q[1],q[0];
rz(-1.91690862178802) q[0];
sx q[0];
rz(1.85146692593629) q[0];
sx q[0];
rz(9.99373058079883) q[0];
rz(0.766674041748047) q[2];
sx q[2];
rz(7.2781230529123) q[2];
sx q[2];
rz(11.9676086664121) q[2];
cx q[2],q[1];
rz(-2.08734488487244) q[1];
sx q[1];
rz(1.99648812611634) q[1];
sx q[1];
rz(11.3991839647214) q[1];
rz(-1.94931447505951) q[3];
sx q[3];
rz(2.64447280962998) q[3];
sx q[3];
rz(11.9395003080289) q[3];
cx q[3],q[2];
rz(0.871716260910034) q[2];
sx q[2];
rz(1.81739524205262) q[2];
sx q[2];
rz(10.3127476930539) q[2];
rz(2.16520833969116) q[3];
sx q[3];
rz(1.72470000584657) q[3];
sx q[3];
rz(11.5606436490934) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.755092144012451) q[0];
sx q[0];
rz(3.21927505930001) q[0];
sx q[0];
rz(9.04715255498096) q[0];
rz(-0.323044210672379) q[1];
sx q[1];
rz(2.47813651164109) q[1];
sx q[1];
rz(10.339948809139) q[1];
cx q[1],q[0];
rz(-0.244190633296967) q[0];
sx q[0];
rz(1.65235975583131) q[0];
sx q[0];
rz(10.672147846214) q[0];
rz(0.56699001789093) q[2];
sx q[2];
rz(5.15123978455598) q[2];
sx q[2];
rz(10.8206287383954) q[2];
cx q[2],q[1];
rz(4.45356369018555) q[1];
sx q[1];
rz(1.88358560402925) q[1];
sx q[1];
rz(9.49696763455077) q[1];
rz(0.876264333724976) q[3];
sx q[3];
rz(1.4834704716974) q[3];
sx q[3];
rz(9.57661830484077) q[3];
cx q[3],q[2];
rz(5.74597549438477) q[2];
sx q[2];
rz(3.30627291102941) q[2];
sx q[2];
rz(5.49266216754123) q[2];
rz(-0.289977133274078) q[3];
sx q[3];
rz(0.732637794809886) q[3];
sx q[3];
rz(13.1837312936704) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.736156046390533) q[0];
sx q[0];
rz(4.20873215992982) q[0];
sx q[0];
rz(9.79181379675075) q[0];
rz(1.23327469825745) q[1];
sx q[1];
rz(5.10925510724122) q[1];
sx q[1];
rz(13.9917621374051) q[1];
cx q[1],q[0];
rz(-0.0504977479577065) q[0];
sx q[0];
rz(4.28005078633363) q[0];
sx q[0];
rz(9.21923022567436) q[0];
rz(-0.0643534660339355) q[2];
sx q[2];
rz(4.26719120343263) q[2];
sx q[2];
rz(10.7264235973279) q[2];
cx q[2],q[1];
rz(-2.51263093948364) q[1];
sx q[1];
rz(5.3169790824228) q[1];
sx q[1];
rz(9.88977090119525) q[1];
rz(2.63787508010864) q[3];
sx q[3];
rz(4.70065382321412) q[3];
sx q[3];
rz(7.06230971812412) q[3];
cx q[3],q[2];
rz(1.16208624839783) q[2];
sx q[2];
rz(4.63381400902803) q[2];
sx q[2];
rz(11.3555955648343) q[2];
rz(-3.13975048065186) q[3];
sx q[3];
rz(3.90709188778932) q[3];
sx q[3];
rz(8.76748511790439) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.84412556886673) q[0];
sx q[0];
rz(3.72589644988114) q[0];
sx q[0];
rz(9.34825971572801) q[0];
rz(0.675297677516937) q[1];
sx q[1];
rz(3.43549820979173) q[1];
sx q[1];
rz(7.03554127215549) q[1];
cx q[1],q[0];
rz(1.99361789226532) q[0];
sx q[0];
rz(2.53208419879014) q[0];
sx q[0];
rz(8.42571077346011) q[0];
rz(0.189578652381897) q[2];
sx q[2];
rz(3.79897329409654) q[2];
sx q[2];
rz(8.3163256406705) q[2];
cx q[2],q[1];
rz(-0.0165377464145422) q[1];
sx q[1];
rz(6.4268948157602) q[1];
sx q[1];
rz(10.0887482523839) q[1];
rz(-2.16032123565674) q[3];
sx q[3];
rz(4.8842889388376) q[3];
sx q[3];
rz(12.9366349935453) q[3];
cx q[3],q[2];
rz(4.07990455627441) q[2];
sx q[2];
rz(4.32180813153321) q[2];
sx q[2];
rz(9.42492445891631) q[2];
rz(1.10953068733215) q[3];
sx q[3];
rz(2.23479816515977) q[3];
sx q[3];
rz(9.72241393326923) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.149348542094231) q[0];
sx q[0];
rz(6.04440966446931) q[0];
sx q[0];
rz(10.4308050632398) q[0];
rz(-0.267931193113327) q[1];
sx q[1];
rz(1.34030989010865) q[1];
sx q[1];
rz(10.7396054029386) q[1];
cx q[1],q[0];
rz(0.490042358636856) q[0];
sx q[0];
rz(2.39391198952729) q[0];
sx q[0];
rz(10.6229210853498) q[0];
rz(-0.133168935775757) q[2];
sx q[2];
rz(1.93179920514161) q[2];
sx q[2];
rz(7.85647222994968) q[2];
cx q[2],q[1];
rz(2.41155743598938) q[1];
sx q[1];
rz(3.95492497284944) q[1];
sx q[1];
rz(8.8398350238721) q[1];
rz(1.57825422286987) q[3];
sx q[3];
rz(2.31416264374787) q[3];
sx q[3];
rz(9.97271445988818) q[3];
cx q[3],q[2];
rz(3.58137059211731) q[2];
sx q[2];
rz(3.68796256382997) q[2];
sx q[2];
rz(9.17930521666213) q[2];
rz(-0.430735111236572) q[3];
sx q[3];
rz(4.23326745827729) q[3];
sx q[3];
rz(13.0436992406766) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.72192871570587) q[0];
sx q[0];
rz(0.992464454966136) q[0];
sx q[0];
rz(9.138096189491) q[0];
rz(0.878967046737671) q[1];
sx q[1];
rz(4.64380136330659) q[1];
sx q[1];
rz(6.67508194445773) q[1];
cx q[1],q[0];
rz(-1.4742728471756) q[0];
sx q[0];
rz(6.84075418313081) q[0];
sx q[0];
rz(9.44049249439641) q[0];
rz(-0.655810177326202) q[2];
sx q[2];
rz(4.82003119786317) q[2];
sx q[2];
rz(12.8937008142392) q[2];
cx q[2],q[1];
rz(-0.380102336406708) q[1];
sx q[1];
rz(5.45156923134858) q[1];
sx q[1];
rz(8.10651323794528) q[1];
rz(2.13042187690735) q[3];
sx q[3];
rz(4.98520997365052) q[3];
sx q[3];
rz(9.61937489210769) q[3];
cx q[3],q[2];
rz(0.6095170378685) q[2];
sx q[2];
rz(4.13459465106065) q[2];
sx q[2];
rz(10.5823230504911) q[2];
rz(-0.550842940807343) q[3];
sx q[3];
rz(4.71939829190309) q[3];
sx q[3];
rz(11.3881640195768) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.206510961055756) q[0];
sx q[0];
rz(5.89932146866853) q[0];
sx q[0];
rz(10.0697213172834) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(3.07462859153748) q[1];
sx q[1];
rz(1.77501288254792) q[1];
sx q[1];
rz(10.412791466705) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(0.930711030960083) q[2];
sx q[2];
rz(0.687340172129222) q[2];
sx q[2];
rz(9.61915648578807) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(2.20636367797852) q[3];
sx q[3];
rz(2.33779081900651) q[3];
sx q[3];
rz(7.04924342631503) q[3];
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