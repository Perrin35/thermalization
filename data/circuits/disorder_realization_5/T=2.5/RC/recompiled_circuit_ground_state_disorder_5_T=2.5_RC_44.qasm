OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.759999) q[0];
sx q[0];
rz(-0.17188369) q[0];
sx q[0];
rz(2.5266393) q[0];
rz(-1.7847269) q[1];
sx q[1];
rz(-1.6109799) q[1];
sx q[1];
rz(-0.15712486) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2295132) q[0];
sx q[0];
rz(-2.3621971) q[0];
sx q[0];
rz(-0.81672658) q[0];
rz(-2.831249) q[2];
sx q[2];
rz(-0.23761286) q[2];
sx q[2];
rz(2.1365191) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8871044) q[1];
sx q[1];
rz(-0.65782065) q[1];
sx q[1];
rz(1.595934) q[1];
x q[2];
rz(-1.8039005) q[3];
sx q[3];
rz(-0.86458428) q[3];
sx q[3];
rz(-1.8854837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5379546) q[2];
sx q[2];
rz(-2.8261638) q[2];
sx q[2];
rz(-1.2856548) q[2];
rz(1.8042608) q[3];
sx q[3];
rz(-2.1888816) q[3];
sx q[3];
rz(0.11346909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86785698) q[0];
sx q[0];
rz(-1.8524167) q[0];
sx q[0];
rz(0.56647545) q[0];
rz(-1.6784338) q[1];
sx q[1];
rz(-0.17371829) q[1];
sx q[1];
rz(-2.4991551) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2447889) q[0];
sx q[0];
rz(-1.519108) q[0];
sx q[0];
rz(1.1079074) q[0];
rz(-1.9011097) q[2];
sx q[2];
rz(-2.1871302) q[2];
sx q[2];
rz(-1.358145) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3088544) q[1];
sx q[1];
rz(-0.843261) q[1];
sx q[1];
rz(3.0738405) q[1];
x q[2];
rz(0.94704721) q[3];
sx q[3];
rz(-0.94175856) q[3];
sx q[3];
rz(-1.8333904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.88175875) q[2];
sx q[2];
rz(-2.3469648) q[2];
sx q[2];
rz(2.3074522) q[2];
rz(0.086890876) q[3];
sx q[3];
rz(-2.0565242) q[3];
sx q[3];
rz(-1.1311919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6777545) q[0];
sx q[0];
rz(-2.0448271) q[0];
sx q[0];
rz(2.0820397) q[0];
rz(-0.73486596) q[1];
sx q[1];
rz(-3.1262472) q[1];
sx q[1];
rz(-2.9954092) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8627166) q[0];
sx q[0];
rz(-1.4116316) q[0];
sx q[0];
rz(1.0907286) q[0];
rz(-pi) q[1];
rz(1.4465815) q[2];
sx q[2];
rz(-1.5746279) q[2];
sx q[2];
rz(1.6506548) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4242658) q[1];
sx q[1];
rz(-1.7317514) q[1];
sx q[1];
rz(2.9907254) q[1];
rz(2.9133818) q[3];
sx q[3];
rz(-1.8442698) q[3];
sx q[3];
rz(0.21450689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.77201468) q[2];
sx q[2];
rz(-2.7507608) q[2];
sx q[2];
rz(2.3005627) q[2];
rz(-0.05989017) q[3];
sx q[3];
rz(-1.6570897) q[3];
sx q[3];
rz(-0.64391518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.197072) q[0];
sx q[0];
rz(-0.50839013) q[0];
sx q[0];
rz(0.12666853) q[0];
rz(1.7238759) q[1];
sx q[1];
rz(-3.1210777) q[1];
sx q[1];
rz(2.1518478) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1331461) q[0];
sx q[0];
rz(-1.5494816) q[0];
sx q[0];
rz(0.7288148) q[0];
x q[1];
rz(1.3700757) q[2];
sx q[2];
rz(-2.7373103) q[2];
sx q[2];
rz(-0.78122369) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9056315) q[1];
sx q[1];
rz(-2.0249019) q[1];
sx q[1];
rz(-2.924861) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7249171) q[3];
sx q[3];
rz(-1.8088578) q[3];
sx q[3];
rz(1.4509089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7030455) q[2];
sx q[2];
rz(-0.35312167) q[2];
sx q[2];
rz(-2.9959196) q[2];
rz(2.6977671) q[3];
sx q[3];
rz(-1.5668818) q[3];
sx q[3];
rz(-1.1183967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1031652) q[0];
sx q[0];
rz(-1.9157836) q[0];
sx q[0];
rz(-0.78381729) q[0];
rz(1.8812284) q[1];
sx q[1];
rz(-3.1385707) q[1];
sx q[1];
rz(0.22617117) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12482551) q[0];
sx q[0];
rz(-1.1065605) q[0];
sx q[0];
rz(-2.8874257) q[0];
x q[1];
rz(-2.1820081) q[2];
sx q[2];
rz(-0.49641616) q[2];
sx q[2];
rz(-0.69545262) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.60701398) q[1];
sx q[1];
rz(-1.4958989) q[1];
sx q[1];
rz(2.043488) q[1];
rz(-pi) q[2];
rz(0.89543476) q[3];
sx q[3];
rz(-0.90762075) q[3];
sx q[3];
rz(1.1010045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1251395) q[2];
sx q[2];
rz(-0.21802248) q[2];
sx q[2];
rz(1.7832635) q[2];
rz(2.6839117) q[3];
sx q[3];
rz(-1.7031368) q[3];
sx q[3];
rz(2.5911736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054475527) q[0];
sx q[0];
rz(-3.1396907) q[0];
sx q[0];
rz(3.0892293) q[0];
rz(-0.26854435) q[1];
sx q[1];
rz(-1.9895376) q[1];
sx q[1];
rz(2.9103738) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5729882) q[0];
sx q[0];
rz(-1.3954961) q[0];
sx q[0];
rz(2.8820287) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8131658) q[2];
sx q[2];
rz(-0.77244416) q[2];
sx q[2];
rz(1.3470115) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0655122) q[1];
sx q[1];
rz(-1.339318) q[1];
sx q[1];
rz(-2.8006018) q[1];
rz(-pi) q[2];
rz(-1.8299425) q[3];
sx q[3];
rz(-2.5185725) q[3];
sx q[3];
rz(-1.3192605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3117567) q[2];
sx q[2];
rz(-2.629403) q[2];
sx q[2];
rz(-2.3066985) q[2];
rz(-3.0541776) q[3];
sx q[3];
rz(-2.0525457) q[3];
sx q[3];
rz(2.1277229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8304623) q[0];
sx q[0];
rz(-0.99263793) q[0];
sx q[0];
rz(-2.2845238) q[0];
rz(-2.3509707) q[1];
sx q[1];
rz(-3.1309083) q[1];
sx q[1];
rz(-1.0766006) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4907376) q[0];
sx q[0];
rz(-1.3044573) q[0];
sx q[0];
rz(-1.2580714) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.78342763) q[2];
sx q[2];
rz(-1.9840709) q[2];
sx q[2];
rz(-2.836712) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8995114) q[1];
sx q[1];
rz(-0.82846009) q[1];
sx q[1];
rz(2.7513585) q[1];
x q[2];
rz(1.7067327) q[3];
sx q[3];
rz(-1.2316717) q[3];
sx q[3];
rz(1.9423758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2662346) q[2];
sx q[2];
rz(-0.42576063) q[2];
sx q[2];
rz(2.2597964) q[2];
rz(-1.734123) q[3];
sx q[3];
rz(-2.2018645) q[3];
sx q[3];
rz(0.56434977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28500104) q[0];
sx q[0];
rz(-0.0017310062) q[0];
sx q[0];
rz(-2.8518387) q[0];
rz(1.2334791) q[1];
sx q[1];
rz(-1.8461485) q[1];
sx q[1];
rz(0.27898702) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.679259) q[0];
sx q[0];
rz(-1.7391608) q[0];
sx q[0];
rz(-3.116045) q[0];
rz(-pi) q[1];
rz(-1.8548429) q[2];
sx q[2];
rz(-1.1835754) q[2];
sx q[2];
rz(-2.7329138) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6304172) q[1];
sx q[1];
rz(-1.8969444) q[1];
sx q[1];
rz(-3.0688968) q[1];
x q[2];
rz(-2.2905228) q[3];
sx q[3];
rz(-1.2382231) q[3];
sx q[3];
rz(0.50902396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4616619) q[2];
sx q[2];
rz(-1.2744023) q[2];
sx q[2];
rz(2.8604841) q[2];
rz(2.5680961) q[3];
sx q[3];
rz(-2.1888013) q[3];
sx q[3];
rz(2.7219685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24116521) q[0];
sx q[0];
rz(-0.92210162) q[0];
sx q[0];
rz(1.4379372) q[0];
rz(-0.18962139) q[1];
sx q[1];
rz(-3.1260243) q[1];
sx q[1];
rz(-1.2356893) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3251117) q[0];
sx q[0];
rz(-2.6345523) q[0];
sx q[0];
rz(-1.965824) q[0];
x q[1];
rz(-0.81053712) q[2];
sx q[2];
rz(-0.97522465) q[2];
sx q[2];
rz(3.0308804) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6744102) q[1];
sx q[1];
rz(-1.7315718) q[1];
sx q[1];
rz(1.2700446) q[1];
rz(-3.1240033) q[3];
sx q[3];
rz(-1.5639728) q[3];
sx q[3];
rz(1.9987035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8427061) q[2];
sx q[2];
rz(-1.7791396) q[2];
sx q[2];
rz(2.5617808) q[2];
rz(1.3696356) q[3];
sx q[3];
rz(-0.42391351) q[3];
sx q[3];
rz(-0.4637318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6746826) q[0];
sx q[0];
rz(-1.5002102) q[0];
sx q[0];
rz(-1.4379733) q[0];
rz(-1.9101608) q[1];
sx q[1];
rz(-2.2302901) q[1];
sx q[1];
rz(-1.4968754) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5176277) q[0];
sx q[0];
rz(-1.8359487) q[0];
sx q[0];
rz(2.318906) q[0];
rz(-pi) q[1];
x q[1];
rz(0.52224285) q[2];
sx q[2];
rz(-1.2127611) q[2];
sx q[2];
rz(-2.4286662) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.10324998) q[1];
sx q[1];
rz(-1.5885067) q[1];
sx q[1];
rz(-1.5679541) q[1];
x q[2];
rz(1.0874463) q[3];
sx q[3];
rz(-1.3030682) q[3];
sx q[3];
rz(-2.8341531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9897495) q[2];
sx q[2];
rz(-1.5785549) q[2];
sx q[2];
rz(2.6452276) q[2];
rz(0.94480354) q[3];
sx q[3];
rz(-0.0050408575) q[3];
sx q[3];
rz(-2.4387824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6481358) q[0];
sx q[0];
rz(-1.428816) q[0];
sx q[0];
rz(-1.1270123) q[0];
rz(1.5588749) q[1];
sx q[1];
rz(-0.3180779) q[1];
sx q[1];
rz(0.19837468) q[1];
rz(1.3140903) q[2];
sx q[2];
rz(-1.5460925) q[2];
sx q[2];
rz(1.8800541) q[2];
rz(-0.54917033) q[3];
sx q[3];
rz(-0.63180379) q[3];
sx q[3];
rz(2.8997811) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
