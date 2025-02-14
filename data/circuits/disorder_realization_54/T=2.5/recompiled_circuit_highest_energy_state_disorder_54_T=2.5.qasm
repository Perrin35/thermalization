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
rz(-0.328313529491425) q[0];
sx q[0];
rz(5.14831224282319) q[0];
sx q[0];
rz(6.84906456469699) q[0];
rz(0.24428753554821) q[1];
sx q[1];
rz(4.92574349244172) q[1];
sx q[1];
rz(9.95903757809802) q[1];
cx q[1],q[0];
rz(-1.78436481952667) q[0];
sx q[0];
rz(4.03816208441789) q[0];
sx q[0];
rz(11.344420528404) q[0];
rz(-5.30731868743896) q[2];
sx q[2];
rz(4.45872870286042) q[2];
sx q[2];
rz(3.81913182734653) q[2];
cx q[2],q[1];
rz(0.457858324050903) q[1];
sx q[1];
rz(3.44301977952058) q[1];
sx q[1];
rz(11.0753910303037) q[1];
rz(2.02796101570129) q[3];
sx q[3];
rz(7.62993827660615) q[3];
sx q[3];
rz(12.2595076322477) q[3];
cx q[3],q[2];
rz(3.12857055664063) q[2];
sx q[2];
rz(3.39209315379197) q[2];
sx q[2];
rz(15.9571475744168) q[2];
rz(4.93919706344604) q[3];
sx q[3];
rz(4.68461433251435) q[3];
sx q[3];
rz(12.6364466905515) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.494211196899414) q[0];
sx q[0];
rz(4.36795380909974) q[0];
sx q[0];
rz(11.5765633344571) q[0];
rz(0.790659427642822) q[1];
sx q[1];
rz(7.05009094079072) q[1];
sx q[1];
rz(9.11312759517833) q[1];
cx q[1],q[0];
rz(2.83833312988281) q[0];
sx q[0];
rz(1.06195345719392) q[0];
sx q[0];
rz(11.8031282186429) q[0];
rz(-4.82059621810913) q[2];
sx q[2];
rz(4.5266993363672) q[2];
sx q[2];
rz(8.64327404498264) q[2];
cx q[2],q[1];
rz(0.0987223014235497) q[1];
sx q[1];
rz(4.47359612782533) q[1];
sx q[1];
rz(7.67678842543765) q[1];
rz(3.75648260116577) q[3];
sx q[3];
rz(1.85281887848909) q[3];
sx q[3];
rz(6.64582250117465) q[3];
cx q[3],q[2];
rz(3.14276838302612) q[2];
sx q[2];
rz(4.30183914502198) q[2];
sx q[2];
rz(10.376161313049) q[2];
rz(-2.913325548172) q[3];
sx q[3];
rz(5.30645862420137) q[3];
sx q[3];
rz(11.2927059888761) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.157628729939461) q[0];
sx q[0];
rz(4.7303263266855) q[0];
sx q[0];
rz(6.1687221288602) q[0];
rz(1.00222563743591) q[1];
sx q[1];
rz(2.71482345660264) q[1];
sx q[1];
rz(12.7282004117887) q[1];
cx q[1],q[0];
rz(-1.96602582931519) q[0];
sx q[0];
rz(5.10095921357209) q[0];
sx q[0];
rz(13.2140407323758) q[0];
rz(4.32349681854248) q[2];
sx q[2];
rz(4.63101735909516) q[2];
sx q[2];
rz(9.97557560204669) q[2];
cx q[2],q[1];
rz(-1.96496498584747) q[1];
sx q[1];
rz(1.33136585553224) q[1];
sx q[1];
rz(8.47670111655399) q[1];
rz(1.97076714038849) q[3];
sx q[3];
rz(0.841639908152171) q[3];
sx q[3];
rz(9.32167344390556) q[3];
cx q[3],q[2];
rz(0.747300446033478) q[2];
sx q[2];
rz(2.69603991706903) q[2];
sx q[2];
rz(6.42394349574252) q[2];
rz(2.58981394767761) q[3];
sx q[3];
rz(5.00915506680543) q[3];
sx q[3];
rz(8.67344877719089) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.93633389472961) q[0];
sx q[0];
rz(6.54922810395295) q[0];
sx q[0];
rz(10.0996693134229) q[0];
rz(-0.154480129480362) q[1];
sx q[1];
rz(4.55063513119752) q[1];
sx q[1];
rz(12.1084048509519) q[1];
cx q[1],q[0];
rz(3.91271591186523) q[0];
sx q[0];
rz(8.24198022683198) q[0];
sx q[0];
rz(14.6674461126248) q[0];
rz(6.88888072967529) q[2];
sx q[2];
rz(6.79155340989167) q[2];
sx q[2];
rz(14.4001917600553) q[2];
cx q[2],q[1];
rz(-1.61246597766876) q[1];
sx q[1];
rz(4.89639559586579) q[1];
sx q[1];
rz(7.47726581095859) q[1];
rz(0.106807552278042) q[3];
sx q[3];
rz(-1.0651009957022) q[3];
sx q[3];
rz(12.4290354013364) q[3];
cx q[3],q[2];
rz(0.133504793047905) q[2];
sx q[2];
rz(3.82183721859986) q[2];
sx q[2];
rz(12.2469234228055) q[2];
rz(1.81146287918091) q[3];
sx q[3];
rz(4.15811673005159) q[3];
sx q[3];
rz(9.98498586415454) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.15956807136536) q[0];
sx q[0];
rz(-4.27790579001372) q[0];
sx q[0];
rz(7.50473318099185) q[0];
rz(-0.232806518673897) q[1];
sx q[1];
rz(2.57229516108567) q[1];
sx q[1];
rz(7.26627371310397) q[1];
cx q[1],q[0];
rz(-0.0966893956065178) q[0];
sx q[0];
rz(4.22466555436189) q[0];
sx q[0];
rz(10.3954483032148) q[0];
rz(2.52091288566589) q[2];
sx q[2];
rz(7.94682231743867) q[2];
sx q[2];
rz(12.862551188461) q[2];
cx q[2],q[1];
rz(-4.14325904846191) q[1];
sx q[1];
rz(7.66772857506806) q[1];
sx q[1];
rz(12.4178275823514) q[1];
rz(3.56391668319702) q[3];
sx q[3];
rz(1.72308007081086) q[3];
sx q[3];
rz(11.4020948171537) q[3];
cx q[3],q[2];
rz(0.382721960544586) q[2];
sx q[2];
rz(7.76754394372041) q[2];
sx q[2];
rz(5.90790102481052) q[2];
rz(2.06598567962646) q[3];
sx q[3];
rz(0.386382492380687) q[3];
sx q[3];
rz(12.0314278364102) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.74068343639374) q[0];
sx q[0];
rz(4.84641161759431) q[0];
sx q[0];
rz(7.98743257521793) q[0];
rz(3.06286931037903) q[1];
sx q[1];
rz(4.51837125618989) q[1];
sx q[1];
rz(12.0182282686154) q[1];
cx q[1],q[0];
rz(2.50544357299805) q[0];
sx q[0];
rz(5.58487311204011) q[0];
sx q[0];
rz(10.6631935596387) q[0];
rz(-0.486058264970779) q[2];
sx q[2];
rz(2.06937387784059) q[2];
sx q[2];
rz(9.28287503718539) q[2];
cx q[2],q[1];
rz(-0.744331479072571) q[1];
sx q[1];
rz(5.10046258767182) q[1];
sx q[1];
rz(11.7089781522672) q[1];
rz(2.32677555084229) q[3];
sx q[3];
rz(5.19235435326631) q[3];
sx q[3];
rz(10.2795093417089) q[3];
cx q[3],q[2];
rz(-1.68462014198303) q[2];
sx q[2];
rz(3.7623012979799) q[2];
sx q[2];
rz(8.74944815634891) q[2];
rz(1.55725491046906) q[3];
sx q[3];
rz(4.9757515509897) q[3];
sx q[3];
rz(11.9492561578672) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.543934404850006) q[0];
sx q[0];
rz(1.18273964722688) q[0];
sx q[0];
rz(9.92077085971042) q[0];
rz(-4.59580326080322) q[1];
sx q[1];
rz(5.22771635850007) q[1];
sx q[1];
rz(8.30801162718936) q[1];
cx q[1],q[0];
rz(1.31632125377655) q[0];
sx q[0];
rz(2.63512352307374) q[0];
sx q[0];
rz(11.783475136749) q[0];
rz(-0.537490427494049) q[2];
sx q[2];
rz(4.47916463215882) q[2];
sx q[2];
rz(7.77213678359195) q[2];
cx q[2],q[1];
rz(-0.361937254667282) q[1];
sx q[1];
rz(1.23142674763734) q[1];
sx q[1];
rz(5.94479486941501) q[1];
rz(1.11478590965271) q[3];
sx q[3];
rz(1.62637618382508) q[3];
sx q[3];
rz(7.35265038012668) q[3];
cx q[3],q[2];
rz(2.3108959197998) q[2];
sx q[2];
rz(2.1572956760698) q[2];
sx q[2];
rz(9.77964872717067) q[2];
rz(0.765346705913544) q[3];
sx q[3];
rz(2.27828750212724) q[3];
sx q[3];
rz(6.37273213862582) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.45277214050293) q[0];
sx q[0];
rz(1.66257623036439) q[0];
sx q[0];
rz(10.1903748869817) q[0];
rz(4.01303291320801) q[1];
sx q[1];
rz(5.51753440697724) q[1];
sx q[1];
rz(8.10203120707675) q[1];
cx q[1],q[0];
rz(0.252155601978302) q[0];
sx q[0];
rz(1.83126345475251) q[0];
sx q[0];
rz(9.3711080647926) q[0];
rz(-8.09041309356689) q[2];
sx q[2];
rz(4.04677501519258) q[2];
sx q[2];
rz(14.7252096891324) q[2];
cx q[2],q[1];
rz(-1.88125276565552) q[1];
sx q[1];
rz(4.69591322739656) q[1];
sx q[1];
rz(4.2040061712186) q[1];
rz(-2.40085816383362) q[3];
sx q[3];
rz(0.301271589594432) q[3];
sx q[3];
rz(11.9190962076108) q[3];
cx q[3],q[2];
rz(-2.48379969596863) q[2];
sx q[2];
rz(4.18435874779756) q[2];
sx q[2];
rz(9.64648973046943) q[2];
rz(1.09295570850372) q[3];
sx q[3];
rz(4.55440548260743) q[3];
sx q[3];
rz(11.8882364988248) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.0738649889826775) q[0];
sx q[0];
rz(4.38818839390809) q[0];
sx q[0];
rz(12.2892102956693) q[0];
rz(-2.39323377609253) q[1];
sx q[1];
rz(-0.447350827855519) q[1];
sx q[1];
rz(11.4230509757917) q[1];
cx q[1],q[0];
rz(1.42440140247345) q[0];
sx q[0];
rz(2.77293080289895) q[0];
sx q[0];
rz(11.9806823491971) q[0];
rz(-2.89509558677673) q[2];
sx q[2];
rz(3.86187491019303) q[2];
sx q[2];
rz(10.9944852352063) q[2];
cx q[2],q[1];
rz(-3.4898955821991) q[1];
sx q[1];
rz(4.86931720574433) q[1];
sx q[1];
rz(6.60422871112033) q[1];
rz(-0.171293452382088) q[3];
sx q[3];
rz(5.06935289700563) q[3];
sx q[3];
rz(10.9937063217084) q[3];
cx q[3],q[2];
rz(-1.91441679000854) q[2];
sx q[2];
rz(4.16514268715913) q[2];
sx q[2];
rz(12.5174035787503) q[2];
rz(-1.06607270240784) q[3];
sx q[3];
rz(5.65480223496492) q[3];
sx q[3];
rz(9.26907379030391) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(3.01860356330872) q[0];
sx q[0];
rz(4.75551572640473) q[0];
sx q[0];
rz(6.26367328166171) q[0];
rz(3.92035841941833) q[1];
sx q[1];
rz(5.59027090867097) q[1];
sx q[1];
rz(10.9876188993375) q[1];
cx q[1],q[0];
rz(3.35854053497314) q[0];
sx q[0];
rz(9.2668835242563) q[0];
sx q[0];
rz(8.82350370883151) q[0];
rz(-0.609930872917175) q[2];
sx q[2];
rz(4.86264458497102) q[2];
sx q[2];
rz(9.40351107380494) q[2];
cx q[2],q[1];
rz(-0.732474684715271) q[1];
sx q[1];
rz(2.15550866921479) q[1];
sx q[1];
rz(11.5225853681485) q[1];
rz(0.164616525173187) q[3];
sx q[3];
rz(0.740629823999949) q[3];
sx q[3];
rz(6.41656420230075) q[3];
cx q[3],q[2];
rz(2.94893002510071) q[2];
sx q[2];
rz(7.14661947091157) q[2];
sx q[2];
rz(15.5424108266751) q[2];
rz(-3.74915838241577) q[3];
sx q[3];
rz(7.56797376473481) q[3];
sx q[3];
rz(8.9564241528432) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.231848686933517) q[0];
sx q[0];
rz(1.59186831315095) q[0];
sx q[0];
rz(11.4770087957303) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(4.57417488098145) q[1];
sx q[1];
rz(3.25353501935536) q[1];
sx q[1];
rz(4.91316125392123) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-0.721257627010345) q[2];
sx q[2];
rz(-0.858697740239553) q[2];
sx q[2];
rz(10.5924998283307) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(0.497292876243591) q[3];
sx q[3];
rz(3.86590632994706) q[3];
sx q[3];
rz(4.15346572398349) q[3];
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
