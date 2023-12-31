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
rz(-0.949375927448273) q[0];
sx q[0];
rz(5.23601427872712) q[0];
sx q[0];
rz(9.49350232481166) q[0];
rz(1.74604952335358) q[1];
sx q[1];
rz(4.67392507393891) q[1];
sx q[1];
rz(8.21640167235538) q[1];
cx q[1],q[0];
rz(2.26573038101196) q[0];
sx q[0];
rz(3.24832809914882) q[0];
sx q[0];
rz(2.69055697917148) q[0];
rz(-1.66628921031952) q[2];
sx q[2];
rz(6.15681591828401) q[2];
sx q[2];
rz(12.1604623556058) q[2];
cx q[2],q[1];
rz(0.613328337669373) q[1];
sx q[1];
rz(4.60168150265748) q[1];
sx q[1];
rz(11.6150078535001) q[1];
rz(-5.44561100006104) q[3];
sx q[3];
rz(2.48237219651277) q[3];
sx q[3];
rz(7.27748725413486) q[3];
cx q[3],q[2];
rz(-1.81571137905121) q[2];
sx q[2];
rz(4.88248923619325) q[2];
sx q[2];
rz(7.95819363593265) q[2];
rz(2.44384336471558) q[3];
sx q[3];
rz(1.10132280190522) q[3];
sx q[3];
rz(10.171938276283) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.117347002029419) q[0];
sx q[0];
rz(4.26924875577027) q[0];
sx q[0];
rz(8.25058422087833) q[0];
rz(-6.4543309211731) q[1];
sx q[1];
rz(1.04479637940461) q[1];
sx q[1];
rz(3.43878886698886) q[1];
cx q[1],q[0];
rz(-3.92431616783142) q[0];
sx q[0];
rz(4.76068571408326) q[0];
sx q[0];
rz(9.50799145399734) q[0];
rz(2.10650944709778) q[2];
sx q[2];
rz(2.3539054115587) q[2];
sx q[2];
rz(10.3860937118451) q[2];
cx q[2],q[1];
rz(-2.61517071723938) q[1];
sx q[1];
rz(2.50896582205827) q[1];
sx q[1];
rz(11.9669804334561) q[1];
rz(-0.808378577232361) q[3];
sx q[3];
rz(4.41794362862641) q[3];
sx q[3];
rz(9.39463677107497) q[3];
cx q[3],q[2];
rz(-0.161956652998924) q[2];
sx q[2];
rz(4.32600704033906) q[2];
sx q[2];
rz(13.1803524255674) q[2];
rz(2.26545143127441) q[3];
sx q[3];
rz(5.61547056038911) q[3];
sx q[3];
rz(10.0618081450383) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.09595942497253) q[0];
sx q[0];
rz(4.99797705014283) q[0];
sx q[0];
rz(9.73241583108112) q[0];
rz(2.39305067062378) q[1];
sx q[1];
rz(3.47314143379266) q[1];
sx q[1];
rz(8.58497062920734) q[1];
cx q[1],q[0];
rz(-2.98999834060669) q[0];
sx q[0];
rz(2.63221058447892) q[0];
sx q[0];
rz(7.24217102526828) q[0];
rz(-1.91018688678741) q[2];
sx q[2];
rz(4.80160680611665) q[2];
sx q[2];
rz(13.3841256856839) q[2];
cx q[2],q[1];
rz(-1.36293435096741) q[1];
sx q[1];
rz(5.16273799737031) q[1];
sx q[1];
rz(6.62309000491306) q[1];
rz(-0.63294130563736) q[3];
sx q[3];
rz(4.85304609139497) q[3];
sx q[3];
rz(7.65141401290103) q[3];
cx q[3],q[2];
rz(2.27411890029907) q[2];
sx q[2];
rz(4.93206623395021) q[2];
sx q[2];
rz(7.60112211703464) q[2];
rz(-1.92582023143768) q[3];
sx q[3];
rz(5.92666712601716) q[3];
sx q[3];
rz(10.9343255519788) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.763862192630768) q[0];
sx q[0];
rz(4.52199176152284) q[0];
sx q[0];
rz(12.1208081006925) q[0];
rz(0.620827794075012) q[1];
sx q[1];
rz(5.03716102440888) q[1];
sx q[1];
rz(11.6007897615354) q[1];
cx q[1],q[0];
rz(1.5774028301239) q[0];
sx q[0];
rz(1.84685424168641) q[0];
sx q[0];
rz(7.07694456576511) q[0];
rz(3.10666155815125) q[2];
sx q[2];
rz(5.55445757706697) q[2];
sx q[2];
rz(6.87446901797458) q[2];
cx q[2],q[1];
rz(7.75356292724609) q[1];
sx q[1];
rz(3.39720025857026) q[1];
sx q[1];
rz(6.40958807467624) q[1];
rz(0.706495404243469) q[3];
sx q[3];
rz(-3.29085764090484) q[3];
sx q[3];
rz(10.0347148537557) q[3];
cx q[3],q[2];
rz(4.46154832839966) q[2];
sx q[2];
rz(3.9039832671457) q[2];
sx q[2];
rz(10.3541033029477) q[2];
rz(0.643541395664215) q[3];
sx q[3];
rz(4.17775181134278) q[3];
sx q[3];
rz(7.28664133547946) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.66995656490326) q[0];
sx q[0];
rz(1.5595373233133) q[0];
sx q[0];
rz(10.7088061332624) q[0];
rz(3.43140268325806) q[1];
sx q[1];
rz(0.739578636484691) q[1];
sx q[1];
rz(13.6145248174588) q[1];
cx q[1],q[0];
rz(-2.78323435783386) q[0];
sx q[0];
rz(4.04951462348039) q[0];
sx q[0];
rz(7.23369739054843) q[0];
rz(-0.628894448280334) q[2];
sx q[2];
rz(3.98558244307572) q[2];
sx q[2];
rz(5.43575165270969) q[2];
cx q[2],q[1];
rz(-0.847612738609314) q[1];
sx q[1];
rz(4.3906456549936) q[1];
sx q[1];
rz(15.8055348157804) q[1];
rz(0.573151230812073) q[3];
sx q[3];
rz(2.37546584208543) q[3];
sx q[3];
rz(11.5123753309171) q[3];
cx q[3],q[2];
rz(5.01228046417236) q[2];
sx q[2];
rz(8.4089619239145) q[2];
sx q[2];
rz(8.52393982409641) q[2];
rz(1.09269309043884) q[3];
sx q[3];
rz(2.13815024693544) q[3];
sx q[3];
rz(13.8004717588346) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.18224048614502) q[0];
sx q[0];
rz(1.93379596074159) q[0];
sx q[0];
rz(10.0208826422612) q[0];
rz(-1.4959362745285) q[1];
sx q[1];
rz(0.897696169214793) q[1];
sx q[1];
rz(11.3214598655622) q[1];
cx q[1],q[0];
rz(5.603600025177) q[0];
sx q[0];
rz(-1.25126585165923) q[0];
sx q[0];
rz(8.89202646016284) q[0];
rz(0.316860526800156) q[2];
sx q[2];
rz(5.17739525635774) q[2];
sx q[2];
rz(13.330176806442) q[2];
cx q[2],q[1];
rz(-0.603312611579895) q[1];
sx q[1];
rz(1.07441750367219) q[1];
sx q[1];
rz(11.7467374563138) q[1];
rz(4.96047878265381) q[3];
sx q[3];
rz(-1.79569420020049) q[3];
sx q[3];
rz(9.5212844222705) q[3];
cx q[3],q[2];
rz(-1.91013765335083) q[2];
sx q[2];
rz(0.698805006342479) q[2];
sx q[2];
rz(9.19981055556937) q[2];
rz(-0.0884300172328949) q[3];
sx q[3];
rz(4.8318499644571) q[3];
sx q[3];
rz(13.0451760053556) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.67696189880371) q[0];
sx q[0];
rz(1.90436521370942) q[0];
sx q[0];
rz(9.14483222960635) q[0];
rz(1.67845582962036) q[1];
sx q[1];
rz(5.01714793046052) q[1];
sx q[1];
rz(9.67746940850421) q[1];
cx q[1],q[0];
rz(1.84259271621704) q[0];
sx q[0];
rz(2.84470662673051) q[0];
sx q[0];
rz(9.58853743075534) q[0];
rz(0.0256645828485489) q[2];
sx q[2];
rz(3.77617588837678) q[2];
sx q[2];
rz(8.88456783293887) q[2];
cx q[2],q[1];
rz(4.74914264678955) q[1];
sx q[1];
rz(4.20708897908265) q[1];
sx q[1];
rz(8.59711948632404) q[1];
rz(2.4239513874054) q[3];
sx q[3];
rz(8.47394910653169) q[3];
sx q[3];
rz(13.0486585855405) q[3];
cx q[3],q[2];
rz(0.320204049348831) q[2];
sx q[2];
rz(5.38951483567292) q[2];
sx q[2];
rz(10.9566490411679) q[2];
rz(1.94847095012665) q[3];
sx q[3];
rz(2.23339978058869) q[3];
sx q[3];
rz(8.33810577391788) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(3.24642610549927) q[0];
sx q[0];
rz(4.90011802514131) q[0];
sx q[0];
rz(10.9108913898389) q[0];
rz(2.87278175354004) q[1];
sx q[1];
rz(5.15322080452973) q[1];
sx q[1];
rz(9.14583126305743) q[1];
cx q[1],q[0];
rz(-4.51418876647949) q[0];
sx q[0];
rz(2.06457522709901) q[0];
sx q[0];
rz(8.70215187071964) q[0];
rz(-0.502119898796082) q[2];
sx q[2];
rz(4.56137481530244) q[2];
sx q[2];
rz(13.0452608823697) q[2];
cx q[2],q[1];
rz(3.76531624794006) q[1];
sx q[1];
rz(4.33807257016236) q[1];
sx q[1];
rz(8.61244455575153) q[1];
rz(-2.56934094429016) q[3];
sx q[3];
rz(4.83748677571351) q[3];
sx q[3];
rz(11.7788893937986) q[3];
cx q[3],q[2];
rz(1.8028609752655) q[2];
sx q[2];
rz(4.60544029076631) q[2];
sx q[2];
rz(7.87585351466342) q[2];
rz(5.23238754272461) q[3];
sx q[3];
rz(7.49014392693574) q[3];
sx q[3];
rz(7.48313460349246) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.67796969413757) q[0];
sx q[0];
rz(4.07517996628816) q[0];
sx q[0];
rz(10.6780004262845) q[0];
rz(0.622500956058502) q[1];
sx q[1];
rz(4.82311180432374) q[1];
sx q[1];
rz(10.5710875749509) q[1];
cx q[1],q[0];
rz(-0.779594004154205) q[0];
sx q[0];
rz(0.75095525582368) q[0];
sx q[0];
rz(7.7509507894437) q[0];
rz(-3.82393431663513) q[2];
sx q[2];
rz(5.25081411202485) q[2];
sx q[2];
rz(10.1741344094197) q[2];
cx q[2],q[1];
rz(-1.51333630084991) q[1];
sx q[1];
rz(5.59207478364045) q[1];
sx q[1];
rz(9.5400610551159) q[1];
rz(-2.20566558837891) q[3];
sx q[3];
rz(1.91626170476014) q[3];
sx q[3];
rz(9.1263754427354) q[3];
cx q[3],q[2];
rz(2.98777318000793) q[2];
sx q[2];
rz(4.17729965050752) q[2];
sx q[2];
rz(10.440609908096) q[2];
rz(0.909795701503754) q[3];
sx q[3];
rz(5.61309471924836) q[3];
sx q[3];
rz(9.05668131112262) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.12904536724091) q[0];
sx q[0];
rz(3.75519517262513) q[0];
sx q[0];
rz(9.41053525115504) q[0];
rz(0.844738900661469) q[1];
sx q[1];
rz(4.09453663428361) q[1];
sx q[1];
rz(11.1847897529523) q[1];
cx q[1],q[0];
rz(-2.25050139427185) q[0];
sx q[0];
rz(1.42824652989442) q[0];
sx q[0];
rz(11.0062569141309) q[0];
rz(2.56270051002502) q[2];
sx q[2];
rz(4.65858724911744) q[2];
sx q[2];
rz(8.66588941811725) q[2];
cx q[2],q[1];
rz(0.0183162931352854) q[1];
sx q[1];
rz(1.01782074769075) q[1];
sx q[1];
rz(10.8265729904096) q[1];
rz(-1.03733456134796) q[3];
sx q[3];
rz(1.76036122639711) q[3];
sx q[3];
rz(10.9520306348722) q[3];
cx q[3],q[2];
rz(1.07667291164398) q[2];
sx q[2];
rz(5.08083966572816) q[2];
sx q[2];
rz(11.5763745069425) q[2];
rz(4.03566360473633) q[3];
sx q[3];
rz(-0.488641349477223) q[3];
sx q[3];
rz(9.65020220576926) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.52790033817291) q[0];
sx q[0];
rz(2.14066651661927) q[0];
sx q[0];
rz(11.3233458757321) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-0.552045226097107) q[1];
sx q[1];
rz(1.79315462906892) q[1];
sx q[1];
rz(10.4102788925092) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(0.865481615066528) q[2];
sx q[2];
rz(4.34242478211457) q[2];
sx q[2];
rz(15.3139967679898) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-1.2726469039917) q[3];
sx q[3];
rz(5.70332876046235) q[3];
sx q[3];
rz(10.8807860374372) q[3];
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
