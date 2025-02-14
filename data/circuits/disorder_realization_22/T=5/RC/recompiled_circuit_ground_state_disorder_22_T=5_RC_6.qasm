OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.30224213) q[0];
sx q[0];
rz(-1.2737162) q[0];
sx q[0];
rz(0.94925517) q[0];
rz(-1.0957837) q[1];
sx q[1];
rz(-0.97691184) q[1];
sx q[1];
rz(1.7916602) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.567141) q[0];
sx q[0];
rz(-1.0691088) q[0];
sx q[0];
rz(-1.6654963) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.737713) q[2];
sx q[2];
rz(-0.35940659) q[2];
sx q[2];
rz(-0.071381005) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7589345) q[1];
sx q[1];
rz(-1.9449502) q[1];
sx q[1];
rz(1.6439245) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1175577) q[3];
sx q[3];
rz(-2.6461678) q[3];
sx q[3];
rz(-1.514243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6076516) q[2];
sx q[2];
rz(-1.093163) q[2];
sx q[2];
rz(-0.1926113) q[2];
rz(-1.8588148) q[3];
sx q[3];
rz(-1.1056113) q[3];
sx q[3];
rz(-2.5530596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.823536) q[0];
sx q[0];
rz(-2.7243491) q[0];
sx q[0];
rz(0.70660025) q[0];
rz(0.63287863) q[1];
sx q[1];
rz(-1.0989847) q[1];
sx q[1];
rz(2.852827) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1624127) q[0];
sx q[0];
rz(-1.5732764) q[0];
sx q[0];
rz(-1.5697877) q[0];
rz(-2.4817746) q[2];
sx q[2];
rz(-2.8809469) q[2];
sx q[2];
rz(-1.6915996) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0367248) q[1];
sx q[1];
rz(-2.782208) q[1];
sx q[1];
rz(-2.4998328) q[1];
x q[2];
rz(0.76008009) q[3];
sx q[3];
rz(-1.2562498) q[3];
sx q[3];
rz(-9.4903367e-05) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.67549813) q[2];
sx q[2];
rz(-2.0530632) q[2];
sx q[2];
rz(-1.4906073) q[2];
rz(-0.7771107) q[3];
sx q[3];
rz(-1.4584352) q[3];
sx q[3];
rz(1.3597663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8260088) q[0];
sx q[0];
rz(-1.6827787) q[0];
sx q[0];
rz(2.7050731) q[0];
rz(-0.79477683) q[1];
sx q[1];
rz(-1.3705285) q[1];
sx q[1];
rz(-0.93903843) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32246958) q[0];
sx q[0];
rz(-1.7717585) q[0];
sx q[0];
rz(-0.51216258) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58231488) q[2];
sx q[2];
rz(-2.2852201) q[2];
sx q[2];
rz(-2.0272875) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3311829) q[1];
sx q[1];
rz(-2.4723834) q[1];
sx q[1];
rz(0.52468467) q[1];
rz(-pi) q[2];
rz(-2.9211798) q[3];
sx q[3];
rz(-2.1142981) q[3];
sx q[3];
rz(1.8956748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7006435) q[2];
sx q[2];
rz(-1.0627397) q[2];
sx q[2];
rz(0.6353333) q[2];
rz(2.3144531) q[3];
sx q[3];
rz(-2.6878036) q[3];
sx q[3];
rz(1.2686096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24469911) q[0];
sx q[0];
rz(-0.06299717) q[0];
sx q[0];
rz(0.52198207) q[0];
rz(2.7105647) q[1];
sx q[1];
rz(-1.5343752) q[1];
sx q[1];
rz(-2.4897051) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70744694) q[0];
sx q[0];
rz(-1.0672671) q[0];
sx q[0];
rz(1.8412526) q[0];
rz(2.2950933) q[2];
sx q[2];
rz(-1.1870459) q[2];
sx q[2];
rz(-0.64696128) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3053599) q[1];
sx q[1];
rz(-1.9499602) q[1];
sx q[1];
rz(-0.01474488) q[1];
x q[2];
rz(-2.1736657) q[3];
sx q[3];
rz(-0.9612007) q[3];
sx q[3];
rz(-0.84679195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7151457) q[2];
sx q[2];
rz(-1.2812252) q[2];
sx q[2];
rz(-0.78488266) q[2];
rz(-2.6134885) q[3];
sx q[3];
rz(-1.2930861) q[3];
sx q[3];
rz(1.2877134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2454979) q[0];
sx q[0];
rz(-2.446785) q[0];
sx q[0];
rz(-0.78952638) q[0];
rz(-0.49742571) q[1];
sx q[1];
rz(-2.1010294) q[1];
sx q[1];
rz(1.8400037) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0987463) q[0];
sx q[0];
rz(-1.7502971) q[0];
sx q[0];
rz(-2.5371206) q[0];
x q[1];
rz(1.5100432) q[2];
sx q[2];
rz(-1.9092602) q[2];
sx q[2];
rz(-2.0450704) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9250946) q[1];
sx q[1];
rz(-1.7431694) q[1];
sx q[1];
rz(2.4045776) q[1];
x q[2];
rz(-0.022707247) q[3];
sx q[3];
rz(-1.2035666) q[3];
sx q[3];
rz(0.87388384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54399458) q[2];
sx q[2];
rz(-1.602729) q[2];
sx q[2];
rz(-2.8653115) q[2];
rz(-0.9497408) q[3];
sx q[3];
rz(-0.26849982) q[3];
sx q[3];
rz(-0.55324078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62040579) q[0];
sx q[0];
rz(-0.036245417) q[0];
sx q[0];
rz(-0.94394839) q[0];
rz(-1.2983407) q[1];
sx q[1];
rz(-1.0669758) q[1];
sx q[1];
rz(-0.7775158) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28431129) q[0];
sx q[0];
rz(-0.67196437) q[0];
sx q[0];
rz(-0.59055968) q[0];
rz(-0.28892679) q[2];
sx q[2];
rz(-2.3863966) q[2];
sx q[2];
rz(-0.64695219) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4247269) q[1];
sx q[1];
rz(-1.9391372) q[1];
sx q[1];
rz(1.4348381) q[1];
rz(-pi) q[2];
x q[2];
rz(1.603807) q[3];
sx q[3];
rz(-1.2668161) q[3];
sx q[3];
rz(-0.98807166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5139318) q[2];
sx q[2];
rz(-1.3769423) q[2];
sx q[2];
rz(-0.39109209) q[2];
rz(-2.5202461) q[3];
sx q[3];
rz(-0.73453271) q[3];
sx q[3];
rz(2.3656316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67052996) q[0];
sx q[0];
rz(-0.10469086) q[0];
sx q[0];
rz(1.712557) q[0];
rz(-2.1757226) q[1];
sx q[1];
rz(-1.8873676) q[1];
sx q[1];
rz(-0.78580725) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98322372) q[0];
sx q[0];
rz(-2.2136723) q[0];
sx q[0];
rz(-1.5063365) q[0];
rz(-pi) q[1];
rz(0.44148032) q[2];
sx q[2];
rz(-0.78462839) q[2];
sx q[2];
rz(-0.048603622) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8829086) q[1];
sx q[1];
rz(-1.6526387) q[1];
sx q[1];
rz(0.46103448) q[1];
rz(-1.2811411) q[3];
sx q[3];
rz(-1.148734) q[3];
sx q[3];
rz(2.6817123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4449731) q[2];
sx q[2];
rz(-0.61903054) q[2];
sx q[2];
rz(0.49717286) q[2];
rz(0.87388006) q[3];
sx q[3];
rz(-2.0495575) q[3];
sx q[3];
rz(1.9397651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1931605) q[0];
sx q[0];
rz(-2.8346297) q[0];
sx q[0];
rz(-2.4440785) q[0];
rz(-2.7613617) q[1];
sx q[1];
rz(-2.0262599) q[1];
sx q[1];
rz(-3.0013705) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2220737) q[0];
sx q[0];
rz(-1.6478879) q[0];
sx q[0];
rz(0.18203966) q[0];
rz(-pi) q[1];
rz(-0.082265286) q[2];
sx q[2];
rz(-2.6994355) q[2];
sx q[2];
rz(1.4210267) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.48988399) q[1];
sx q[1];
rz(-1.9700642) q[1];
sx q[1];
rz(0.7266161) q[1];
rz(-1.7626552) q[3];
sx q[3];
rz(-1.975112) q[3];
sx q[3];
rz(-1.9003001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2908638) q[2];
sx q[2];
rz(-1.9242492) q[2];
sx q[2];
rz(0.49120894) q[2];
rz(-2.9344007) q[3];
sx q[3];
rz(-0.62316337) q[3];
sx q[3];
rz(1.4108968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17289138) q[0];
sx q[0];
rz(-1.7270813) q[0];
sx q[0];
rz(-0.15039314) q[0];
rz(-2.4405759) q[1];
sx q[1];
rz(-2.0471408) q[1];
sx q[1];
rz(1.4775803) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1883066) q[0];
sx q[0];
rz(-0.64670783) q[0];
sx q[0];
rz(0.20215277) q[0];
rz(-pi) q[1];
rz(-1.4544746) q[2];
sx q[2];
rz(-1.7309628) q[2];
sx q[2];
rz(0.016906658) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7733998) q[1];
sx q[1];
rz(-2.9559694) q[1];
sx q[1];
rz(2.5297861) q[1];
rz(-pi) q[2];
rz(1.2712237) q[3];
sx q[3];
rz(-0.75787395) q[3];
sx q[3];
rz(1.8488499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.56665862) q[2];
sx q[2];
rz(-0.40859544) q[2];
sx q[2];
rz(-1.5236141) q[2];
rz(0.59897113) q[3];
sx q[3];
rz(-2.6409769) q[3];
sx q[3];
rz(-2.9536501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69342518) q[0];
sx q[0];
rz(-0.42891112) q[0];
sx q[0];
rz(-1.5018139) q[0];
rz(-0.93694726) q[1];
sx q[1];
rz(-1.0277156) q[1];
sx q[1];
rz(1.017259) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27148025) q[0];
sx q[0];
rz(-1.9471629) q[0];
sx q[0];
rz(0.41618213) q[0];
x q[1];
rz(-1.4849365) q[2];
sx q[2];
rz(-1.6756264) q[2];
sx q[2];
rz(1.3121999) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5084475) q[1];
sx q[1];
rz(-0.26403759) q[1];
sx q[1];
rz(-1.4435084) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3689092) q[3];
sx q[3];
rz(-1.1491778) q[3];
sx q[3];
rz(-3.0516171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0246058) q[2];
sx q[2];
rz(-2.459343) q[2];
sx q[2];
rz(-1.5218706) q[2];
rz(-1.4874602) q[3];
sx q[3];
rz(-1.6493075) q[3];
sx q[3];
rz(-1.3967995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7617154) q[0];
sx q[0];
rz(-1.7424255) q[0];
sx q[0];
rz(-2.4095834) q[0];
rz(-1.4749745) q[1];
sx q[1];
rz(-1.2154308) q[1];
sx q[1];
rz(-2.348127) q[1];
rz(-0.39948612) q[2];
sx q[2];
rz(-2.1151092) q[2];
sx q[2];
rz(-1.5110037) q[2];
rz(-0.31994896) q[3];
sx q[3];
rz(-1.2816164) q[3];
sx q[3];
rz(1.5641227) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
