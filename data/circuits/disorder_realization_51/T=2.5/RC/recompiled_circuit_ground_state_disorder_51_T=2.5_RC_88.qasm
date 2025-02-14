OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.92276031) q[0];
sx q[0];
rz(3.2050686) q[0];
sx q[0];
rz(10.795074) q[0];
rz(0.81075794) q[1];
sx q[1];
rz(-1.4501362) q[1];
sx q[1];
rz(-2.9210747) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9618586) q[0];
sx q[0];
rz(-2.4027236) q[0];
sx q[0];
rz(1.4703077) q[0];
rz(-pi) q[1];
rz(-2.0283255) q[2];
sx q[2];
rz(-1.2731247) q[2];
sx q[2];
rz(3.0106017) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.045956417) q[1];
sx q[1];
rz(-1.7929309) q[1];
sx q[1];
rz(-1.5820811) q[1];
rz(-pi) q[2];
rz(1.0706717) q[3];
sx q[3];
rz(-2.7410738) q[3];
sx q[3];
rz(2.5409215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.084879547) q[2];
sx q[2];
rz(-3.1352545) q[2];
sx q[2];
rz(-2.32178) q[2];
rz(0.020717185) q[3];
sx q[3];
rz(-2.8262704) q[3];
sx q[3];
rz(-1.0516385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3725975) q[0];
sx q[0];
rz(-0.1854493) q[0];
sx q[0];
rz(-0.73082596) q[0];
rz(-1.7164879) q[1];
sx q[1];
rz(-2.4162879) q[1];
sx q[1];
rz(-0.46754974) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.088652164) q[0];
sx q[0];
rz(-3.0346617) q[0];
sx q[0];
rz(1.8631975) q[0];
x q[1];
rz(3.0655083) q[2];
sx q[2];
rz(-2.7307163) q[2];
sx q[2];
rz(-2.8166688) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7035571) q[1];
sx q[1];
rz(-2.272637) q[1];
sx q[1];
rz(1.7662394) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2456513) q[3];
sx q[3];
rz(-2.3811901) q[3];
sx q[3];
rz(1.5261306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9759489) q[2];
sx q[2];
rz(-0.011198137) q[2];
sx q[2];
rz(0.33640081) q[2];
rz(-2.8339556) q[3];
sx q[3];
rz(-3.1308789) q[3];
sx q[3];
rz(0.78905869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12863185) q[0];
sx q[0];
rz(-2.9427981) q[0];
sx q[0];
rz(3.0096753) q[0];
rz(1.4384653) q[1];
sx q[1];
rz(-1.4142282) q[1];
sx q[1];
rz(-1.4836813) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8298873) q[0];
sx q[0];
rz(-2.1714806) q[0];
sx q[0];
rz(-1.2313103) q[0];
rz(1.5882115) q[2];
sx q[2];
rz(-1.5952186) q[2];
sx q[2];
rz(-1.2877854) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10569948) q[1];
sx q[1];
rz(-1.4669682) q[1];
sx q[1];
rz(3.0894214) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68861945) q[3];
sx q[3];
rz(-3.0303749) q[3];
sx q[3];
rz(0.86513774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5710473) q[2];
sx q[2];
rz(-1.8071625) q[2];
sx q[2];
rz(-1.6837616) q[2];
rz(2.3633862) q[3];
sx q[3];
rz(-3.1359378) q[3];
sx q[3];
rz(-0.2125423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0197765) q[0];
sx q[0];
rz(-1.7541405) q[0];
sx q[0];
rz(2.8507267) q[0];
rz(1.5485171) q[1];
sx q[1];
rz(-2.3389108) q[1];
sx q[1];
rz(3.1150418) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7598068) q[0];
sx q[0];
rz(-1.2273754) q[0];
sx q[0];
rz(1.9476101) q[0];
x q[1];
rz(1.6241811) q[2];
sx q[2];
rz(-2.0211271) q[2];
sx q[2];
rz(0.93102294) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.091152126) q[1];
sx q[1];
rz(-0.9602957) q[1];
sx q[1];
rz(0.97368413) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7351355) q[3];
sx q[3];
rz(-3.1222635) q[3];
sx q[3];
rz(3.1301182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.87963858) q[2];
sx q[2];
rz(-3.1231572) q[2];
sx q[2];
rz(-0.43799841) q[2];
rz(-2.0336464) q[3];
sx q[3];
rz(-0.017939311) q[3];
sx q[3];
rz(-1.3560791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9172025) q[0];
sx q[0];
rz(-0.49773911) q[0];
sx q[0];
rz(2.9955731) q[0];
rz(3.0085425) q[1];
sx q[1];
rz(-2.0818043) q[1];
sx q[1];
rz(2.6859247) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9216292) q[0];
sx q[0];
rz(-1.5908003) q[0];
sx q[0];
rz(-0.09394333) q[0];
rz(-pi) q[1];
rz(-2.5495569) q[2];
sx q[2];
rz(-1.8285654) q[2];
sx q[2];
rz(-1.3883615) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.86381493) q[1];
sx q[1];
rz(-1.3425273) q[1];
sx q[1];
rz(2.9991908) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5862053) q[3];
sx q[3];
rz(-1.3432627) q[3];
sx q[3];
rz(-1.7209382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2666152) q[2];
sx q[2];
rz(-0.019860331) q[2];
sx q[2];
rz(1.88545) q[2];
rz(-0.52136326) q[3];
sx q[3];
rz(-2.8837236) q[3];
sx q[3];
rz(2.7969587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94619036) q[0];
sx q[0];
rz(-2.7552216) q[0];
sx q[0];
rz(-0.56889164) q[0];
rz(-0.16828123) q[1];
sx q[1];
rz(-1.5852837) q[1];
sx q[1];
rz(0.010244244) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83547384) q[0];
sx q[0];
rz(-0.17806192) q[0];
sx q[0];
rz(-2.8639069) q[0];
rz(3.1222759) q[2];
sx q[2];
rz(-1.3903119) q[2];
sx q[2];
rz(1.3767124) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0371984) q[1];
sx q[1];
rz(-1.5362829) q[1];
sx q[1];
rz(-1.6134279) q[1];
rz(0.041268392) q[3];
sx q[3];
rz(-2.5070243) q[3];
sx q[3];
rz(0.41555017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.16113082) q[2];
sx q[2];
rz(-0.22393301) q[2];
sx q[2];
rz(1.0978318) q[2];
rz(0.3499507) q[3];
sx q[3];
rz(-1.2772468) q[3];
sx q[3];
rz(-2.2866975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2646645) q[0];
sx q[0];
rz(-2.9783037) q[0];
sx q[0];
rz(0.40633416) q[0];
rz(-1.3111275) q[1];
sx q[1];
rz(-0.51478148) q[1];
sx q[1];
rz(-3.0313671) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7007164) q[0];
sx q[0];
rz(-1.5982107) q[0];
sx q[0];
rz(3.0423881) q[0];
x q[1];
rz(1.5761779) q[2];
sx q[2];
rz(-1.5670243) q[2];
sx q[2];
rz(1.9413858) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4496042) q[1];
sx q[1];
rz(-1.9983728) q[1];
sx q[1];
rz(3.1319992) q[1];
rz(-pi) q[2];
rz(-0.87894999) q[3];
sx q[3];
rz(-0.96262299) q[3];
sx q[3];
rz(-1.6339621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.280764) q[2];
sx q[2];
rz(-0.0088366652) q[2];
sx q[2];
rz(0.78738085) q[2];
rz(0.27960882) q[3];
sx q[3];
rz(-3.0968554) q[3];
sx q[3];
rz(-2.4217822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1396007) q[0];
sx q[0];
rz(-0.75374341) q[0];
sx q[0];
rz(-1.7107704) q[0];
rz(0.17499533) q[1];
sx q[1];
rz(-0.7376968) q[1];
sx q[1];
rz(1.7134679) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.437182) q[0];
sx q[0];
rz(-1.5488273) q[0];
sx q[0];
rz(-3.1385413) q[0];
x q[1];
rz(-2.0753209) q[2];
sx q[2];
rz(-3.1331535) q[2];
sx q[2];
rz(1.5167459) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0261532) q[1];
sx q[1];
rz(-2.8255531) q[1];
sx q[1];
rz(-1.0277961) q[1];
x q[2];
rz(-2.8599127) q[3];
sx q[3];
rz(-1.3550645) q[3];
sx q[3];
rz(-0.40116596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.78508198) q[2];
sx q[2];
rz(-2.3729237) q[2];
sx q[2];
rz(-1.53995) q[2];
rz(-0.29940638) q[3];
sx q[3];
rz(-1.6573903) q[3];
sx q[3];
rz(-0.33777344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47141075) q[0];
sx q[0];
rz(-3.1306559) q[0];
sx q[0];
rz(2.6543044) q[0];
rz(1.4393073) q[1];
sx q[1];
rz(-0.7936365) q[1];
sx q[1];
rz(0.38584858) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7283994) q[0];
sx q[0];
rz(-1.5110713) q[0];
sx q[0];
rz(-2.7483536) q[0];
rz(-3.0294777) q[2];
sx q[2];
rz(-1.9939853) q[2];
sx q[2];
rz(-1.8433169) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6035406) q[1];
sx q[1];
rz(-1.6038938) q[1];
sx q[1];
rz(2.4500928) q[1];
x q[2];
rz(1.5351546) q[3];
sx q[3];
rz(-1.9058133) q[3];
sx q[3];
rz(-0.50594508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7172598) q[2];
sx q[2];
rz(-0.029742664) q[2];
sx q[2];
rz(-2.7359803) q[2];
rz(-0.82617104) q[3];
sx q[3];
rz(-3.0982389) q[3];
sx q[3];
rz(-0.9894754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9721603) q[0];
sx q[0];
rz(-2.5763474) q[0];
sx q[0];
rz(1.7461079) q[0];
rz(-1.6204429) q[1];
sx q[1];
rz(-2.4902159) q[1];
sx q[1];
rz(1.7924538) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17332213) q[0];
sx q[0];
rz(-1.7815457) q[0];
sx q[0];
rz(0.83619946) q[0];
rz(-pi) q[1];
rz(-1.8113891) q[2];
sx q[2];
rz(-1.4052011) q[2];
sx q[2];
rz(2.9179433) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6031462) q[1];
sx q[1];
rz(-1.6308404) q[1];
sx q[1];
rz(0.014343099) q[1];
rz(-pi) q[2];
rz(-0.46197148) q[3];
sx q[3];
rz(-3.1079227) q[3];
sx q[3];
rz(-1.8505423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.32538432) q[2];
sx q[2];
rz(-3.1367229) q[2];
sx q[2];
rz(-1.7083141) q[2];
rz(1.2895182) q[3];
sx q[3];
rz(-3.0472445) q[3];
sx q[3];
rz(-1.2963699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0636487) q[0];
sx q[0];
rz(-2.1049121) q[0];
sx q[0];
rz(-2.5459469) q[0];
rz(0.044364914) q[1];
sx q[1];
rz(-0.27635074) q[1];
sx q[1];
rz(1.9377294) q[1];
rz(-0.2169256) q[2];
sx q[2];
rz(-1.1033162) q[2];
sx q[2];
rz(0.14324506) q[2];
rz(-1.9808922) q[3];
sx q[3];
rz(-1.6312508) q[3];
sx q[3];
rz(-1.4283258) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
