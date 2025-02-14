OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.15314161777496) q[0];
sx q[0];
rz(2.32601312001283) q[0];
sx q[0];
rz(10.1829671025197) q[0];
rz(-0.0125013329088688) q[1];
sx q[1];
rz(4.97954359849031) q[1];
sx q[1];
rz(10.9951675891797) q[1];
cx q[1],q[0];
rz(0.859382927417755) q[0];
sx q[0];
rz(3.87891683180863) q[0];
sx q[0];
rz(10.4468365669171) q[0];
rz(1.63816368579865) q[2];
sx q[2];
rz(4.31359568436677) q[2];
sx q[2];
rz(10.1138066410939) q[2];
cx q[2],q[1];
rz(0.0581522881984711) q[1];
sx q[1];
rz(4.65526893933351) q[1];
sx q[1];
rz(9.78562784790202) q[1];
rz(0.309239029884338) q[3];
sx q[3];
rz(4.65494898160035) q[3];
sx q[3];
rz(9.74784693717166) q[3];
cx q[3],q[2];
rz(1.61151731014252) q[2];
sx q[2];
rz(3.13341834780807) q[2];
sx q[2];
rz(8.98407891987964) q[2];
rz(0.0797444507479668) q[3];
sx q[3];
rz(3.14169714155799) q[3];
sx q[3];
rz(10.548925971977) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.175286203622818) q[0];
sx q[0];
rz(3.19845506002242) q[0];
sx q[0];
rz(9.25677598117992) q[0];
rz(0.0202782303094864) q[1];
sx q[1];
rz(3.44985756476457) q[1];
sx q[1];
rz(10.9613717555921) q[1];
cx q[1],q[0];
rz(-0.288233011960983) q[0];
sx q[0];
rz(3.46368873317773) q[0];
sx q[0];
rz(10.249180173866) q[0];
rz(-0.241705507040024) q[2];
sx q[2];
rz(3.1520134789222) q[2];
sx q[2];
rz(8.07041738032504) q[2];
cx q[2],q[1];
rz(0.550243377685547) q[1];
sx q[1];
rz(3.1463920605653) q[1];
sx q[1];
rz(9.65582502483531) q[1];
rz(1.36130034923553) q[3];
sx q[3];
rz(3.09835463215644) q[3];
sx q[3];
rz(10.3033429145734) q[3];
cx q[3],q[2];
rz(0.345404714345932) q[2];
sx q[2];
rz(4.05350104172761) q[2];
sx q[2];
rz(10.8324421405713) q[2];
rz(1.04508543014526) q[3];
sx q[3];
rz(3.09206913610036) q[3];
sx q[3];
rz(9.15277234315082) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.831809103488922) q[0];
sx q[0];
rz(4.11826929648454) q[0];
sx q[0];
rz(9.9858265876691) q[0];
rz(0.276841193437576) q[1];
sx q[1];
rz(3.12871530571068) q[1];
sx q[1];
rz(10.7325867175977) q[1];
cx q[1],q[0];
rz(0.553717434406281) q[0];
sx q[0];
rz(3.44128382404382) q[0];
sx q[0];
rz(9.28553951381847) q[0];
rz(1.57146632671356) q[2];
sx q[2];
rz(3.14898193779448) q[2];
sx q[2];
rz(10.4875808715741) q[2];
cx q[2],q[1];
rz(1.77371156215668) q[1];
sx q[1];
rz(3.71786806185777) q[1];
sx q[1];
rz(9.31726602315112) q[1];
rz(-0.550753116607666) q[3];
sx q[3];
rz(3.8982003053003) q[3];
sx q[3];
rz(10.1570173859517) q[3];
cx q[3],q[2];
rz(1.34485590457916) q[2];
sx q[2];
rz(3.14171075064736) q[2];
sx q[2];
rz(10.0141618609349) q[2];
rz(1.89925897121429) q[3];
sx q[3];
rz(3.15396038920666) q[3];
sx q[3];
rz(10.7850047111432) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.00684833433479071) q[0];
sx q[0];
rz(3.65341046650941) q[0];
sx q[0];
rz(11.2177108287732) q[0];
rz(-0.00665612472221255) q[1];
sx q[1];
rz(4.46278062661225) q[1];
sx q[1];
rz(9.45827877371713) q[1];
cx q[1],q[0];
rz(1.63146996498108) q[0];
sx q[0];
rz(4.06646052201326) q[0];
sx q[0];
rz(10.7786062717359) q[0];
rz(-0.00449677323922515) q[2];
sx q[2];
rz(4.59275403817231) q[2];
sx q[2];
rz(9.8436015009801) q[2];
cx q[2],q[1];
rz(-0.120431452989578) q[1];
sx q[1];
rz(3.41068527300889) q[1];
sx q[1];
rz(9.38653758390948) q[1];
rz(0.636568188667297) q[3];
sx q[3];
rz(3.34498815436895) q[3];
sx q[3];
rz(10.4682864904325) q[3];
cx q[3],q[2];
rz(0.0324307046830654) q[2];
sx q[2];
rz(3.13506317150826) q[2];
sx q[2];
rz(9.72342652677699) q[2];
rz(1.2700891494751) q[3];
sx q[3];
rz(3.12543694314594) q[3];
sx q[3];
rz(9.37158616109892) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.77151596546173) q[0];
sx q[0];
rz(4.6560367663675) q[0];
sx q[0];
rz(9.87162763475581) q[0];
rz(-0.186137855052948) q[1];
sx q[1];
rz(3.2033998948806) q[1];
sx q[1];
rz(11.1532326698224) q[1];
cx q[1],q[0];
rz(-0.227677747607231) q[0];
sx q[0];
rz(3.37488150795037) q[0];
sx q[0];
rz(9.27880336939498) q[0];
rz(-0.422225564718246) q[2];
sx q[2];
rz(1.76864257653291) q[2];
sx q[2];
rz(10.0374488592069) q[2];
cx q[2],q[1];
rz(0.912921011447906) q[1];
sx q[1];
rz(3.20861442585523) q[1];
sx q[1];
rz(10.1633317232053) q[1];
rz(0.666543006896973) q[3];
sx q[3];
rz(2.78911408980424) q[3];
sx q[3];
rz(10.1721019506375) q[3];
cx q[3],q[2];
rz(-0.83429217338562) q[2];
sx q[2];
rz(4.69363072712953) q[2];
sx q[2];
rz(8.92609748839542) q[2];
rz(0.577916085720062) q[3];
sx q[3];
rz(2.65798649390275) q[3];
sx q[3];
rz(10.0214839935224) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.961082875728607) q[0];
sx q[0];
rz(4.26235762436921) q[0];
sx q[0];
rz(9.7711910366933) q[0];
rz(-0.601804256439209) q[1];
sx q[1];
rz(4.70249107678468) q[1];
sx q[1];
rz(10.1782817602079) q[1];
cx q[1],q[0];
rz(1.84180557727814) q[0];
sx q[0];
rz(3.40046841104562) q[0];
sx q[0];
rz(10.2635650396268) q[0];
rz(0.591624021530151) q[2];
sx q[2];
rz(3.2756925543123) q[2];
sx q[2];
rz(9.45023037529691) q[2];
cx q[2],q[1];
rz(0.063976377248764) q[1];
sx q[1];
rz(4.29815569718415) q[1];
sx q[1];
rz(10.1300497412603) q[1];
rz(0.512291848659515) q[3];
sx q[3];
rz(3.32817593415315) q[3];
sx q[3];
rz(9.45761824994489) q[3];
cx q[3],q[2];
rz(0.570091307163239) q[2];
sx q[2];
rz(3.13810205774429) q[2];
sx q[2];
rz(11.0422558545987) q[2];
rz(0.120247185230255) q[3];
sx q[3];
rz(3.13829385547946) q[3];
sx q[3];
rz(9.96252385377094) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.268712490797043) q[0];
sx q[0];
rz(4.06554332573945) q[0];
sx q[0];
rz(9.20351687668964) q[0];
rz(1.46067535877228) q[1];
sx q[1];
rz(4.07497725089128) q[1];
sx q[1];
rz(9.50207872538968) q[1];
cx q[1],q[0];
rz(0.54228264093399) q[0];
sx q[0];
rz(3.099838269996) q[0];
sx q[0];
rz(9.38970890491411) q[0];
rz(0.494281470775604) q[2];
sx q[2];
rz(3.15167255898053) q[2];
sx q[2];
rz(9.14220098256274) q[2];
cx q[2],q[1];
rz(-0.128800228238106) q[1];
sx q[1];
rz(2.95230527420575) q[1];
sx q[1];
rz(10.6066409110944) q[1];
rz(-0.213996186852455) q[3];
sx q[3];
rz(4.83711019356782) q[3];
sx q[3];
rz(9.00492954849407) q[3];
cx q[3],q[2];
rz(1.34954249858856) q[2];
sx q[2];
rz(3.1304065511846) q[2];
sx q[2];
rz(10.3847388982694) q[2];
rz(0.329444259405136) q[3];
sx q[3];
rz(3.1335148123377) q[3];
sx q[3];
rz(10.2750649213712) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.946211695671082) q[0];
sx q[0];
rz(3.7546952684694) q[0];
sx q[0];
rz(9.53405929207011) q[0];
rz(-0.373361349105835) q[1];
sx q[1];
rz(3.9513135274225) q[1];
sx q[1];
rz(11.3359088659207) q[1];
cx q[1],q[0];
rz(0.627235352993011) q[0];
sx q[0];
rz(4.11472568114335) q[0];
sx q[0];
rz(10.2282648444097) q[0];
rz(1.37353789806366) q[2];
sx q[2];
rz(4.52567354043061) q[2];
sx q[2];
rz(9.44183157793387) q[2];
cx q[2],q[1];
rz(0.444303929805756) q[1];
sx q[1];
rz(3.07378008415038) q[1];
sx q[1];
rz(9.25528184174701) q[1];
rz(0.6508749127388) q[3];
sx q[3];
rz(4.9616711457544) q[3];
sx q[3];
rz(10.5951398372571) q[3];
cx q[3],q[2];
rz(1.56651484966278) q[2];
sx q[2];
rz(4.37735501130159) q[2];
sx q[2];
rz(10.7533727645795) q[2];
rz(1.39684200286865) q[3];
sx q[3];
rz(3.14533305157023) q[3];
sx q[3];
rz(10.4466571569364) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.0638655722141266) q[0];
sx q[0];
rz(4.58461049397523) q[0];
sx q[0];
rz(9.99778691529437) q[0];
rz(-0.308144629001617) q[1];
sx q[1];
rz(3.55146351655061) q[1];
sx q[1];
rz(11.5589749574582) q[1];
cx q[1],q[0];
rz(0.389781683683395) q[0];
sx q[0];
rz(3.03490188916261) q[0];
sx q[0];
rz(9.47707152589365) q[0];
rz(-0.307290852069855) q[2];
sx q[2];
rz(3.2858319153362) q[2];
sx q[2];
rz(9.76950923203632) q[2];
cx q[2],q[1];
rz(-0.633576214313507) q[1];
sx q[1];
rz(3.23417707731063) q[1];
sx q[1];
rz(13.0212688207547) q[1];
rz(0.419024348258972) q[3];
sx q[3];
rz(3.10858266626532) q[3];
sx q[3];
rz(9.21075359582111) q[3];
cx q[3],q[2];
rz(-1.31743001937866) q[2];
sx q[2];
rz(3.76835933526094) q[2];
sx q[2];
rz(9.81473583578273) q[2];
rz(0.0690848752856255) q[3];
sx q[3];
rz(3.13248129946227) q[3];
sx q[3];
rz(9.75631222724124) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.0201417151838541) q[0];
sx q[0];
rz(3.89262518485124) q[0];
sx q[0];
rz(9.9107191324155) q[0];
rz(-0.871569752693176) q[1];
sx q[1];
rz(4.97531715233857) q[1];
sx q[1];
rz(10.9190466165464) q[1];
cx q[1],q[0];
rz(0.450099736452103) q[0];
sx q[0];
rz(2.1039823611551) q[0];
sx q[0];
rz(9.78109342455074) q[0];
rz(-0.0441877953708172) q[2];
sx q[2];
rz(2.18515971501405) q[2];
sx q[2];
rz(9.4938440605919) q[2];
cx q[2],q[1];
rz(-0.680771291255951) q[1];
sx q[1];
rz(3.5938288291269) q[1];
sx q[1];
rz(10.2534379720609) q[1];
rz(-0.0419765971601009) q[3];
sx q[3];
rz(1.68383613427217) q[3];
sx q[3];
rz(9.93247113227054) q[3];
cx q[3],q[2];
rz(1.57470524311066) q[2];
sx q[2];
rz(3.18437538121874) q[2];
sx q[2];
rz(9.39289956390067) q[2];
rz(0.778306424617767) q[3];
sx q[3];
rz(3.14840820816393) q[3];
sx q[3];
rz(9.72028086184665) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.30059146881104) q[0];
sx q[0];
rz(3.38837394316728) q[0];
sx q[0];
rz(9.2672640889804) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(1.48386406898499) q[1];
sx q[1];
rz(1.62247064908082) q[1];
sx q[1];
rz(10.7620982885282) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(1.81277930736542) q[2];
sx q[2];
rz(4.74063602288301) q[2];
sx q[2];
rz(10.8586631774823) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(0.361503332853317) q[3];
sx q[3];
rz(5.07187953789765) q[3];
sx q[3];
rz(10.1544627308766) q[3];
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
