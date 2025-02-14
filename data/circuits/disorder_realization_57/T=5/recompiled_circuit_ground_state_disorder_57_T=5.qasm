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
rz(-0.398555546998978) q[0];
sx q[0];
rz(2.27502128680284) q[0];
sx q[0];
rz(9.08781430720493) q[0];
rz(0.269518673419952) q[1];
sx q[1];
rz(4.08428469498689) q[1];
sx q[1];
rz(10.6842912197034) q[1];
cx q[1],q[0];
rz(-1.24526059627533) q[0];
sx q[0];
rz(4.17515710194642) q[0];
sx q[0];
rz(11.8974840402524) q[0];
rz(-0.210026547312737) q[2];
sx q[2];
rz(3.98632571299607) q[2];
sx q[2];
rz(12.4619695901792) q[2];
cx q[2],q[1];
rz(-0.0776905417442322) q[1];
sx q[1];
rz(3.69614848692948) q[1];
sx q[1];
rz(11.1853066444318) q[1];
rz(-3.74101948738098) q[3];
sx q[3];
rz(4.36482623417909) q[3];
sx q[3];
rz(11.7336427926938) q[3];
cx q[3],q[2];
rz(0.172446146607399) q[2];
sx q[2];
rz(4.70598021348054) q[2];
sx q[2];
rz(10.658609366409) q[2];
rz(1.5305757522583) q[3];
sx q[3];
rz(1.73500350316102) q[3];
sx q[3];
rz(9.98988053797885) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.81561404466629) q[0];
sx q[0];
rz(4.10666474898393) q[0];
sx q[0];
rz(9.66876766680881) q[0];
rz(-1.98323929309845) q[1];
sx q[1];
rz(5.26901760895784) q[1];
sx q[1];
rz(8.97643843888446) q[1];
cx q[1],q[0];
rz(0.876544237136841) q[0];
sx q[0];
rz(3.66565397580201) q[0];
sx q[0];
rz(8.21249339579746) q[0];
rz(2.50162839889526) q[2];
sx q[2];
rz(3.99242863257463) q[2];
sx q[2];
rz(8.9065122961919) q[2];
cx q[2],q[1];
rz(4.06544017791748) q[1];
sx q[1];
rz(3.8780706842714) q[1];
sx q[1];
rz(5.58174059390231) q[1];
rz(0.710426032543182) q[3];
sx q[3];
rz(5.60521641572053) q[3];
sx q[3];
rz(10.6600044727246) q[3];
cx q[3],q[2];
rz(1.19251942634583) q[2];
sx q[2];
rz(3.10188803647692) q[2];
sx q[2];
rz(11.7555370092313) q[2];
rz(0.00909463223069906) q[3];
sx q[3];
rz(4.75706413586671) q[3];
sx q[3];
rz(9.10724419950649) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.41191649436951) q[0];
sx q[0];
rz(4.93878999550874) q[0];
sx q[0];
rz(11.6186311006467) q[0];
rz(2.24322271347046) q[1];
sx q[1];
rz(3.85463789303834) q[1];
sx q[1];
rz(9.63112645446464) q[1];
cx q[1],q[0];
rz(1.9131098985672) q[0];
sx q[0];
rz(3.89646539290483) q[0];
sx q[0];
rz(9.41362160294458) q[0];
rz(-0.401385426521301) q[2];
sx q[2];
rz(3.96094438632066) q[2];
sx q[2];
rz(8.13027582167789) q[2];
cx q[2],q[1];
rz(0.628132998943329) q[1];
sx q[1];
rz(5.16467157204682) q[1];
sx q[1];
rz(10.7663221120755) q[1];
rz(-0.0469303354620934) q[3];
sx q[3];
rz(5.03841820557649) q[3];
sx q[3];
rz(9.66707748769923) q[3];
cx q[3],q[2];
rz(-0.729489386081696) q[2];
sx q[2];
rz(2.14958426554734) q[2];
sx q[2];
rz(11.5821308851163) q[2];
rz(-0.497368305921555) q[3];
sx q[3];
rz(4.40091112454469) q[3];
sx q[3];
rz(10.1697793960492) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.36181485652924) q[0];
sx q[0];
rz(3.23346959252889) q[0];
sx q[0];
rz(9.63111490606471) q[0];
rz(-1.67747175693512) q[1];
sx q[1];
rz(0.762814911203929) q[1];
sx q[1];
rz(8.79864243268176) q[1];
cx q[1],q[0];
rz(1.04476058483124) q[0];
sx q[0];
rz(4.01486847003038) q[0];
sx q[0];
rz(7.38927433489963) q[0];
rz(0.00176499981898814) q[2];
sx q[2];
rz(4.8151754458719) q[2];
sx q[2];
rz(10.1202542543332) q[2];
cx q[2],q[1];
rz(1.5765677690506) q[1];
sx q[1];
rz(3.89301642973954) q[1];
sx q[1];
rz(8.59394416808292) q[1];
rz(1.47186470031738) q[3];
sx q[3];
rz(5.09739187558229) q[3];
sx q[3];
rz(10.510307288162) q[3];
cx q[3],q[2];
rz(0.765544056892395) q[2];
sx q[2];
rz(3.91930601199205) q[2];
sx q[2];
rz(10.4145366907041) q[2];
rz(1.03422129154205) q[3];
sx q[3];
rz(5.11415484746034) q[3];
sx q[3];
rz(8.89921048878833) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.70222103595734) q[0];
sx q[0];
rz(4.86344233353669) q[0];
sx q[0];
rz(11.4035361766736) q[0];
rz(0.167126640677452) q[1];
sx q[1];
rz(4.66685250599916) q[1];
sx q[1];
rz(8.74940983056232) q[1];
cx q[1],q[0];
rz(0.576083123683929) q[0];
sx q[0];
rz(4.49414161046083) q[0];
sx q[0];
rz(10.3254602909009) q[0];
rz(-0.379470646381378) q[2];
sx q[2];
rz(4.39539781411225) q[2];
sx q[2];
rz(7.73068950175449) q[2];
cx q[2],q[1];
rz(0.310842216014862) q[1];
sx q[1];
rz(4.28435138066346) q[1];
sx q[1];
rz(10.6145064592282) q[1];
rz(-0.393463164567947) q[3];
sx q[3];
rz(4*pi/5) q[3];
sx q[3];
rz(10.2846148967664) q[3];
cx q[3],q[2];
rz(-1.21657884120941) q[2];
sx q[2];
rz(4.86068967183168) q[2];
sx q[2];
rz(12.0588125944059) q[2];
rz(-0.0192858595401049) q[3];
sx q[3];
rz(3.9365960081392) q[3];
sx q[3];
rz(8.20544395445987) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.602922022342682) q[0];
sx q[0];
rz(2.51410731871659) q[0];
sx q[0];
rz(8.72311625479862) q[0];
rz(0.642746686935425) q[1];
sx q[1];
rz(4.30093625386293) q[1];
sx q[1];
rz(10.7458154916684) q[1];
cx q[1],q[0];
rz(-0.821339905261993) q[0];
sx q[0];
rz(2.51393154461915) q[0];
sx q[0];
rz(9.44629401750072) q[0];
rz(1.21395802497864) q[2];
sx q[2];
rz(5.38044133980805) q[2];
sx q[2];
rz(10.2673843264501) q[2];
cx q[2],q[1];
rz(-0.630326449871063) q[1];
sx q[1];
rz(4.18512669404084) q[1];
sx q[1];
rz(11.4265331983487) q[1];
rz(-1.85943603515625) q[3];
sx q[3];
rz(2.30967900355393) q[3];
sx q[3];
rz(12.1120021104734) q[3];
cx q[3],q[2];
rz(2.84644365310669) q[2];
sx q[2];
rz(3.62290603120858) q[2];
sx q[2];
rz(8.78499904870197) q[2];
rz(2.48447370529175) q[3];
sx q[3];
rz(4.63403967221315) q[3];
sx q[3];
rz(9.1415580868642) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.0397286340594292) q[0];
sx q[0];
rz(4.42777732213075) q[0];
sx q[0];
rz(10.2052379012029) q[0];
rz(1.23774933815002) q[1];
sx q[1];
rz(1.84338084061677) q[1];
sx q[1];
rz(9.26073808073207) q[1];
cx q[1],q[0];
rz(0.177241608500481) q[0];
sx q[0];
rz(3.72105661233003) q[0];
sx q[0];
rz(9.86816162466213) q[0];
rz(0.443497061729431) q[2];
sx q[2];
rz(1.94645431836183) q[2];
sx q[2];
rz(10.7083226203839) q[2];
cx q[2],q[1];
rz(-0.611728131771088) q[1];
sx q[1];
rz(4.19698146184022) q[1];
sx q[1];
rz(11.5881995916288) q[1];
rz(0.458254247903824) q[3];
sx q[3];
rz(2.30831048090989) q[3];
sx q[3];
rz(9.02089629172488) q[3];
cx q[3],q[2];
rz(-0.268745452165604) q[2];
sx q[2];
rz(5.08612492878968) q[2];
sx q[2];
rz(8.54409656523868) q[2];
rz(-1.66704797744751) q[3];
sx q[3];
rz(2.07379523118074) q[3];
sx q[3];
rz(11.1517724752347) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.473776519298553) q[0];
sx q[0];
rz(3.33033901651437) q[0];
sx q[0];
rz(11.4371759653012) q[0];
rz(0.953436255455017) q[1];
sx q[1];
rz(2.2013180573755) q[1];
sx q[1];
rz(10.0133014678876) q[1];
cx q[1],q[0];
rz(-0.557109653949738) q[0];
sx q[0];
rz(3.77721652586991) q[0];
sx q[0];
rz(8.90560302733585) q[0];
rz(-0.332576274871826) q[2];
sx q[2];
rz(1.78754487832124) q[2];
sx q[2];
rz(9.76931399702235) q[2];
cx q[2],q[1];
rz(-0.578916311264038) q[1];
sx q[1];
rz(4.06956348021562) q[1];
sx q[1];
rz(9.82670865058109) q[1];
rz(2.91300821304321) q[3];
sx q[3];
rz(2.91160634358461) q[3];
sx q[3];
rz(9.40543218179747) q[3];
cx q[3],q[2];
rz(-0.15527655184269) q[2];
sx q[2];
rz(4.14940384228761) q[2];
sx q[2];
rz(11.4741611242215) q[2];
rz(0.808825969696045) q[3];
sx q[3];
rz(4.22328880627687) q[3];
sx q[3];
rz(7.72735998629733) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.587619304656982) q[0];
sx q[0];
rz(2.06461361249024) q[0];
sx q[0];
rz(8.75314191579028) q[0];
rz(-0.492975860834122) q[1];
sx q[1];
rz(4.02235076029832) q[1];
sx q[1];
rz(10.6643255710523) q[1];
cx q[1],q[0];
rz(-0.171044707298279) q[0];
sx q[0];
rz(4.47515490849549) q[0];
sx q[0];
rz(10.9982397317807) q[0];
rz(-0.00995306391268969) q[2];
sx q[2];
rz(4.51516524155671) q[2];
sx q[2];
rz(8.27338061331912) q[2];
cx q[2],q[1];
rz(1.08604919910431) q[1];
sx q[1];
rz(2.14308390219743) q[1];
sx q[1];
rz(9.94889811276599) q[1];
rz(-1.40268528461456) q[3];
sx q[3];
rz(3.91765764554078) q[3];
sx q[3];
rz(10.123962378494) q[3];
cx q[3],q[2];
rz(-0.129418313503265) q[2];
sx q[2];
rz(4.61538270314271) q[2];
sx q[2];
rz(9.97134248017474) q[2];
rz(2.21873450279236) q[3];
sx q[3];
rz(5.6542216857248) q[3];
sx q[3];
rz(11.4198468685071) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.657611429691315) q[0];
sx q[0];
rz(3.60168627102906) q[0];
sx q[0];
rz(9.85381031631633) q[0];
rz(1.28267621994019) q[1];
sx q[1];
rz(4.50692096550996) q[1];
sx q[1];
rz(9.50603691338702) q[1];
cx q[1],q[0];
rz(-0.0704393982887268) q[0];
sx q[0];
rz(2.95500463445718) q[0];
sx q[0];
rz(10.5576131105344) q[0];
rz(-1.60333752632141) q[2];
sx q[2];
rz(4.25250449975068) q[2];
sx q[2];
rz(8.64877698420688) q[2];
cx q[2],q[1];
rz(0.729662001132965) q[1];
sx q[1];
rz(4.03747245867784) q[1];
sx q[1];
rz(10.4647487163465) q[1];
rz(1.10544514656067) q[3];
sx q[3];
rz(2.02922597725923) q[3];
sx q[3];
rz(7.398718571655) q[3];
cx q[3],q[2];
rz(-0.73605865240097) q[2];
sx q[2];
rz(4.07249769766862) q[2];
sx q[2];
rz(11.4547137975614) q[2];
rz(-1.97236549854279) q[3];
sx q[3];
rz(3.85477581818635) q[3];
sx q[3];
rz(9.56776035427257) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.812035083770752) q[0];
sx q[0];
rz(5.1700977404886) q[0];
sx q[0];
rz(11.434858775131) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(0.685521364212036) q[1];
sx q[1];
rz(4.75027945836122) q[1];
sx q[1];
rz(9.68185586332485) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-0.248726382851601) q[2];
sx q[2];
rz(4.46041539509828) q[2];
sx q[2];
rz(11.0471431970517) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(0.694139301776886) q[3];
sx q[3];
rz(1.92928078969056) q[3];
sx q[3];
rz(11.5478052854459) q[3];
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
