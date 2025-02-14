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
rz(-3.0781167) q[0];
sx q[0];
rz(1.3702962) q[0];
rz(0.81075794) q[1];
sx q[1];
rz(-1.4501362) q[1];
sx q[1];
rz(-2.9210747) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.044174319) q[0];
sx q[0];
rz(-0.83651354) q[0];
sx q[0];
rz(-0.091139779) q[0];
rz(-pi) q[1];
rz(-2.8120997) q[2];
sx q[2];
rz(-1.1348083) q[2];
sx q[2];
rz(1.5832251) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.097135492) q[1];
sx q[1];
rz(-2.9191764) q[1];
sx q[1];
rz(-3.0916721) q[1];
x q[2];
rz(2.070921) q[3];
sx q[3];
rz(-0.4005188) q[3];
sx q[3];
rz(-0.60067117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.084879547) q[2];
sx q[2];
rz(-3.1352545) q[2];
sx q[2];
rz(2.32178) q[2];
rz(-3.1208755) q[3];
sx q[3];
rz(-2.8262704) q[3];
sx q[3];
rz(2.0899541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3725975) q[0];
sx q[0];
rz(-0.1854493) q[0];
sx q[0];
rz(0.73082596) q[0];
rz(1.4251047) q[1];
sx q[1];
rz(-0.72530472) q[1];
sx q[1];
rz(-2.6740429) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38263812) q[0];
sx q[0];
rz(-1.6731723) q[0];
sx q[0];
rz(3.1106614) q[0];
rz(-pi) q[1];
rz(-0.076084332) q[2];
sx q[2];
rz(-2.7307163) q[2];
sx q[2];
rz(-2.8166688) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.135957) q[1];
sx q[1];
rz(-1.7196481) q[1];
sx q[1];
rz(-0.71133228) q[1];
rz(-1.2456513) q[3];
sx q[3];
rz(-0.76040254) q[3];
sx q[3];
rz(1.5261306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9759489) q[2];
sx q[2];
rz(-3.1303945) q[2];
sx q[2];
rz(-2.8051918) q[2];
rz(-2.8339556) q[3];
sx q[3];
rz(-0.010713723) q[3];
sx q[3];
rz(2.352534) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0129608) q[0];
sx q[0];
rz(-0.19879453) q[0];
sx q[0];
rz(-0.13191731) q[0];
rz(1.4384653) q[1];
sx q[1];
rz(-1.7273644) q[1];
sx q[1];
rz(1.4836813) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0795201) q[0];
sx q[0];
rz(-1.8490845) q[0];
sx q[0];
rz(0.62836439) q[0];
x q[1];
rz(-1.5882115) q[2];
sx q[2];
rz(-1.546374) q[2];
sx q[2];
rz(-1.2877854) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4705088) q[1];
sx q[1];
rz(-1.5189063) q[1];
sx q[1];
rz(1.6747649) q[1];
rz(0.68861945) q[3];
sx q[3];
rz(-0.11121777) q[3];
sx q[3];
rz(2.2764549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5710473) q[2];
sx q[2];
rz(-1.8071625) q[2];
sx q[2];
rz(-1.6837616) q[2];
rz(0.7782065) q[3];
sx q[3];
rz(-0.0056548803) q[3];
sx q[3];
rz(2.9290504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1218162) q[0];
sx q[0];
rz(-1.3874522) q[0];
sx q[0];
rz(0.29086599) q[0];
rz(-1.5485171) q[1];
sx q[1];
rz(-0.80268186) q[1];
sx q[1];
rz(3.1150418) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3817859) q[0];
sx q[0];
rz(-1.9142173) q[0];
sx q[0];
rz(-1.1939826) q[0];
rz(0.10992421) q[2];
sx q[2];
rz(-0.4532686) q[2];
sx q[2];
rz(-2.088415) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77912731) q[1];
sx q[1];
rz(-0.82632768) q[1];
sx q[1];
rz(2.4646548) q[1];
rz(-pi) q[2];
rz(1.5517275) q[3];
sx q[3];
rz(-1.5676343) q[3];
sx q[3];
rz(-1.3950128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2619541) q[2];
sx q[2];
rz(-3.1231572) q[2];
sx q[2];
rz(2.7035942) q[2];
rz(2.0336464) q[3];
sx q[3];
rz(-0.017939311) q[3];
sx q[3];
rz(1.3560791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22439013) q[0];
sx q[0];
rz(-0.49773911) q[0];
sx q[0];
rz(2.9955731) q[0];
rz(-0.13305013) q[1];
sx q[1];
rz(-2.0818043) q[1];
sx q[1];
rz(-0.45566794) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3489482) q[0];
sx q[0];
rz(-1.6647208) q[0];
sx q[0];
rz(1.5908888) q[0];
x q[1];
rz(-2.5495569) q[2];
sx q[2];
rz(-1.8285654) q[2];
sx q[2];
rz(-1.3883615) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2777777) q[1];
sx q[1];
rz(-1.3425273) q[1];
sx q[1];
rz(0.14240188) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0751441) q[3];
sx q[3];
rz(-2.9135468) q[3];
sx q[3];
rz(-1.4888637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8749775) q[2];
sx q[2];
rz(-0.019860331) q[2];
sx q[2];
rz(1.2561426) q[2];
rz(0.52136326) q[3];
sx q[3];
rz(-0.25786906) q[3];
sx q[3];
rz(-0.34463394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94619036) q[0];
sx q[0];
rz(-0.38637105) q[0];
sx q[0];
rz(-2.572701) q[0];
rz(2.9733114) q[1];
sx q[1];
rz(-1.5852837) q[1];
sx q[1];
rz(-3.1313484) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1327577) q[0];
sx q[0];
rz(-1.6193701) q[0];
sx q[0];
rz(0.17137613) q[0];
x q[1];
rz(-3.1222759) q[2];
sx q[2];
rz(-1.7512808) q[2];
sx q[2];
rz(1.3767124) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.10439428) q[1];
sx q[1];
rz(-1.5362829) q[1];
sx q[1];
rz(-1.5281648) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63416173) q[3];
sx q[3];
rz(-1.5463357) q[3];
sx q[3];
rz(2.0195877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.16113082) q[2];
sx q[2];
rz(-2.9176596) q[2];
sx q[2];
rz(1.0978318) q[2];
rz(-2.791642) q[3];
sx q[3];
rz(-1.2772468) q[3];
sx q[3];
rz(-2.2866975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8769281) q[0];
sx q[0];
rz(-2.9783037) q[0];
sx q[0];
rz(0.40633416) q[0];
rz(-1.3111275) q[1];
sx q[1];
rz(-0.51478148) q[1];
sx q[1];
rz(-3.0313671) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12719181) q[0];
sx q[0];
rz(-1.6699635) q[0];
sx q[0];
rz(1.5983461) q[0];
rz(-pi) q[1];
rz(-2.1821457) q[2];
sx q[2];
rz(-0.0065718647) q[2];
sx q[2];
rz(2.159664) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4264723) q[1];
sx q[1];
rz(-2.7139152) q[1];
sx q[1];
rz(-1.5497472) q[1];
x q[2];
rz(-0.87894999) q[3];
sx q[3];
rz(-0.96262299) q[3];
sx q[3];
rz(1.5076306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8608287) q[2];
sx q[2];
rz(-0.0088366652) q[2];
sx q[2];
rz(-2.3542118) q[2];
rz(0.27960882) q[3];
sx q[3];
rz(-3.0968554) q[3];
sx q[3];
rz(0.71981049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0019919458) q[0];
sx q[0];
rz(-2.3878492) q[0];
sx q[0];
rz(1.7107704) q[0];
rz(2.9665973) q[1];
sx q[1];
rz(-0.7376968) q[1];
sx q[1];
rz(1.4281248) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0080453) q[0];
sx q[0];
rz(-1.5677457) q[0];
sx q[0];
rz(1.5488272) q[0];
x q[1];
rz(1.0662717) q[2];
sx q[2];
rz(-3.1331535) q[2];
sx q[2];
rz(1.5167459) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0261532) q[1];
sx q[1];
rz(-0.31603957) q[1];
sx q[1];
rz(2.1137966) q[1];
x q[2];
rz(-1.7950906) q[3];
sx q[3];
rz(-1.845775) q[3];
sx q[3];
rz(2.0338273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.78508198) q[2];
sx q[2];
rz(-0.76866895) q[2];
sx q[2];
rz(1.6016426) q[2];
rz(0.29940638) q[3];
sx q[3];
rz(-1.4842024) q[3];
sx q[3];
rz(2.8038192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6701819) q[0];
sx q[0];
rz(-3.1306559) q[0];
sx q[0];
rz(0.48728824) q[0];
rz(-1.4393073) q[1];
sx q[1];
rz(-0.7936365) q[1];
sx q[1];
rz(2.7557441) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13284616) q[0];
sx q[0];
rz(-1.1782968) q[0];
sx q[0];
rz(-1.6354435) q[0];
rz(-3.0294777) q[2];
sx q[2];
rz(-1.1476074) q[2];
sx q[2];
rz(1.8433169) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5380521) q[1];
sx q[1];
rz(-1.5376989) q[1];
sx q[1];
rz(0.69149986) q[1];
x q[2];
rz(-3.039592) q[3];
sx q[3];
rz(-0.33683646) q[3];
sx q[3];
rz(-2.5276195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4243329) q[2];
sx q[2];
rz(-3.11185) q[2];
sx q[2];
rz(-0.40561238) q[2];
rz(0.82617104) q[3];
sx q[3];
rz(-0.043353733) q[3];
sx q[3];
rz(2.1521173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9721603) q[0];
sx q[0];
rz(-2.5763474) q[0];
sx q[0];
rz(1.3954847) q[0];
rz(-1.5211498) q[1];
sx q[1];
rz(-0.65137678) q[1];
sx q[1];
rz(-1.3491389) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9717108) q[0];
sx q[0];
rz(-0.75877178) q[0];
sx q[0];
rz(-1.8797329) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95913943) q[2];
sx q[2];
rz(-0.29116071) q[2];
sx q[2];
rz(0.75555246) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7730798) q[1];
sx q[1];
rz(-0.061731438) q[1];
sx q[1];
rz(-1.8049989) q[1];
rz(2.6796212) q[3];
sx q[3];
rz(-3.1079227) q[3];
sx q[3];
rz(1.2910503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8162083) q[2];
sx q[2];
rz(-3.1367229) q[2];
sx q[2];
rz(1.4332786) q[2];
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
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0636487) q[0];
sx q[0];
rz(-1.0366806) q[0];
sx q[0];
rz(0.59564577) q[0];
rz(-3.0972277) q[1];
sx q[1];
rz(-0.27635074) q[1];
sx q[1];
rz(1.9377294) q[1];
rz(-1.9738214) q[2];
sx q[2];
rz(-2.6296305) q[2];
sx q[2];
rz(0.59811022) q[2];
rz(-1.4201319) q[3];
sx q[3];
rz(-0.41427783) q[3];
sx q[3];
rz(-2.8610341) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
