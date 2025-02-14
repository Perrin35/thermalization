OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0889283) q[0];
sx q[0];
rz(-1.1298236) q[0];
sx q[0];
rz(3.1273754) q[0];
rz(0.05834236) q[1];
sx q[1];
rz(-2.3925233) q[1];
sx q[1];
rz(0.60325375) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5644419) q[0];
sx q[0];
rz(-1.9214726) q[0];
sx q[0];
rz(-0.68564424) q[0];
rz(-pi) q[1];
x q[1];
rz(1.649734) q[2];
sx q[2];
rz(-2.4853155) q[2];
sx q[2];
rz(-2.5205534) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3989714) q[1];
sx q[1];
rz(-1.4150029) q[1];
sx q[1];
rz(0.69424082) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.18666191) q[3];
sx q[3];
rz(-0.87282729) q[3];
sx q[3];
rz(1.1503134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3732036) q[2];
sx q[2];
rz(-1.5378121) q[2];
sx q[2];
rz(-1.0962037) q[2];
rz(-1.0783892) q[3];
sx q[3];
rz(-2.6069141) q[3];
sx q[3];
rz(-0.53511846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.815149) q[0];
sx q[0];
rz(-2.09477) q[0];
sx q[0];
rz(-0.4775508) q[0];
rz(2.1304456) q[1];
sx q[1];
rz(-2.5944581) q[1];
sx q[1];
rz(-1.0844885) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80974977) q[0];
sx q[0];
rz(-1.6879775) q[0];
sx q[0];
rz(-1.940889) q[0];
x q[1];
rz(-1.162503) q[2];
sx q[2];
rz(-1.5973685) q[2];
sx q[2];
rz(-1.0607189) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1275) q[1];
sx q[1];
rz(-2.1154101) q[1];
sx q[1];
rz(-1.3633534) q[1];
rz(2.4020477) q[3];
sx q[3];
rz(-1.7266183) q[3];
sx q[3];
rz(2.9062944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.61887211) q[2];
sx q[2];
rz(-2.7174157) q[2];
sx q[2];
rz(2.0514533) q[2];
rz(0.56385309) q[3];
sx q[3];
rz(-0.68958759) q[3];
sx q[3];
rz(3.1088945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76145935) q[0];
sx q[0];
rz(-0.021012336) q[0];
sx q[0];
rz(1.6047961) q[0];
rz(-1.8236632) q[1];
sx q[1];
rz(-1.0937966) q[1];
sx q[1];
rz(2.7440548) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5814272) q[0];
sx q[0];
rz(-2.4545963) q[0];
sx q[0];
rz(1.2233748) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0084559) q[2];
sx q[2];
rz(-2.5175885) q[2];
sx q[2];
rz(-0.74988264) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.51008114) q[1];
sx q[1];
rz(-0.95247522) q[1];
sx q[1];
rz(0.31429283) q[1];
rz(-2.0769172) q[3];
sx q[3];
rz(-1.1017513) q[3];
sx q[3];
rz(-2.7739547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.90618187) q[2];
sx q[2];
rz(-1.1778888) q[2];
sx q[2];
rz(0.57514352) q[2];
rz(0.67342526) q[3];
sx q[3];
rz(-1.6005102) q[3];
sx q[3];
rz(-1.5967691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6509318) q[0];
sx q[0];
rz(-2.0199825) q[0];
sx q[0];
rz(-0.90890539) q[0];
rz(-0.017223651) q[1];
sx q[1];
rz(-2.6135018) q[1];
sx q[1];
rz(-0.012103279) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4151013) q[0];
sx q[0];
rz(-2.4612777) q[0];
sx q[0];
rz(0.52796396) q[0];
x q[1];
rz(2.1075756) q[2];
sx q[2];
rz(-2.1488948) q[2];
sx q[2];
rz(2.5930717) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3084532) q[1];
sx q[1];
rz(-2.1670682) q[1];
sx q[1];
rz(0.27383974) q[1];
x q[2];
rz(1.6161259) q[3];
sx q[3];
rz(-1.6974276) q[3];
sx q[3];
rz(1.8243197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.072731344) q[2];
sx q[2];
rz(-2.8673745) q[2];
sx q[2];
rz(-0.38632986) q[2];
rz(2.7522411) q[3];
sx q[3];
rz(-1.6081622) q[3];
sx q[3];
rz(-2.1336335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35578457) q[0];
sx q[0];
rz(-2.4171827) q[0];
sx q[0];
rz(-0.32387787) q[0];
rz(-0.38662275) q[1];
sx q[1];
rz(-1.7522248) q[1];
sx q[1];
rz(1.5868384) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51832685) q[0];
sx q[0];
rz(-0.73438209) q[0];
sx q[0];
rz(0.096303864) q[0];
rz(-pi) q[1];
rz(-1.0350758) q[2];
sx q[2];
rz(-1.8958099) q[2];
sx q[2];
rz(-2.6485788) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.388098) q[1];
sx q[1];
rz(-0.83067229) q[1];
sx q[1];
rz(0.64708935) q[1];
rz(-pi) q[2];
rz(0.79894508) q[3];
sx q[3];
rz(-1.9513592) q[3];
sx q[3];
rz(2.9137672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3891478) q[2];
sx q[2];
rz(-2.2989595) q[2];
sx q[2];
rz(3.0774934) q[2];
rz(2.5373503) q[3];
sx q[3];
rz(-1.7490381) q[3];
sx q[3];
rz(1.909168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9391249) q[0];
sx q[0];
rz(-2.5305643) q[0];
sx q[0];
rz(-2.2623999) q[0];
rz(-2.7912256) q[1];
sx q[1];
rz(-1.8354974) q[1];
sx q[1];
rz(2.9894357) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50987303) q[0];
sx q[0];
rz(-1.3001967) q[0];
sx q[0];
rz(2.4218049) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9486676) q[2];
sx q[2];
rz(-2.0461651) q[2];
sx q[2];
rz(-1.0032636) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.85206807) q[1];
sx q[1];
rz(-1.7729323) q[1];
sx q[1];
rz(-2.7006671) q[1];
x q[2];
rz(1.9655439) q[3];
sx q[3];
rz(-1.9874032) q[3];
sx q[3];
rz(-0.28232251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1019613) q[2];
sx q[2];
rz(-0.32932082) q[2];
sx q[2];
rz(-2.9015818) q[2];
rz(0.65252423) q[3];
sx q[3];
rz(-1.1381166) q[3];
sx q[3];
rz(1.2751381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0535102) q[0];
sx q[0];
rz(-1.2367915) q[0];
sx q[0];
rz(-2.2834593) q[0];
rz(1.3872604) q[1];
sx q[1];
rz(-0.45448449) q[1];
sx q[1];
rz(-0.33285704) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8008566) q[0];
sx q[0];
rz(-1.8854257) q[0];
sx q[0];
rz(-0.97140177) q[0];
x q[1];
rz(0.24999745) q[2];
sx q[2];
rz(-1.6277024) q[2];
sx q[2];
rz(-0.075486334) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.27204) q[1];
sx q[1];
rz(-1.7764047) q[1];
sx q[1];
rz(-0.42821347) q[1];
rz(-0.44705963) q[3];
sx q[3];
rz(-2.4363228) q[3];
sx q[3];
rz(1.7763607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3420928) q[2];
sx q[2];
rz(-0.19793333) q[2];
sx q[2];
rz(-2.0905154) q[2];
rz(-1.4970695) q[3];
sx q[3];
rz(-1.1727419) q[3];
sx q[3];
rz(2.3433949) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.368211) q[0];
sx q[0];
rz(-2.1147275) q[0];
sx q[0];
rz(2.2494702) q[0];
rz(2.6719773) q[1];
sx q[1];
rz(-0.62925595) q[1];
sx q[1];
rz(2.3687252) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5443916) q[0];
sx q[0];
rz(-2.3344451) q[0];
sx q[0];
rz(-0.76274422) q[0];
x q[1];
rz(2.8373884) q[2];
sx q[2];
rz(-0.72663621) q[2];
sx q[2];
rz(-1.8223423) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.23456) q[1];
sx q[1];
rz(-1.2508243) q[1];
sx q[1];
rz(2.6522308) q[1];
rz(-1.8065213) q[3];
sx q[3];
rz(-1.8060214) q[3];
sx q[3];
rz(-2.1947088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.61432394) q[2];
sx q[2];
rz(-1.9169151) q[2];
sx q[2];
rz(0.69250715) q[2];
rz(-0.45037371) q[3];
sx q[3];
rz(-1.2053442) q[3];
sx q[3];
rz(1.6374121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5699128) q[0];
sx q[0];
rz(-2.3501861) q[0];
sx q[0];
rz(-0.94863844) q[0];
rz(2.0750849) q[1];
sx q[1];
rz(-0.98943168) q[1];
sx q[1];
rz(-1.8705286) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97590896) q[0];
sx q[0];
rz(-1.6635487) q[0];
sx q[0];
rz(1.3646056) q[0];
x q[1];
rz(-1.5367212) q[2];
sx q[2];
rz(-1.5928895) q[2];
sx q[2];
rz(1.3349229) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8066602) q[1];
sx q[1];
rz(-2.1927823) q[1];
sx q[1];
rz(-2.9750664) q[1];
rz(-pi) q[2];
rz(-0.44235559) q[3];
sx q[3];
rz(-1.6887553) q[3];
sx q[3];
rz(0.30260146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5991685) q[2];
sx q[2];
rz(-1.8343265) q[2];
sx q[2];
rz(3.001281) q[2];
rz(0.03242759) q[3];
sx q[3];
rz(-0.92493886) q[3];
sx q[3];
rz(-1.9353346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62823826) q[0];
sx q[0];
rz(-1.6980549) q[0];
sx q[0];
rz(-2.3640609) q[0];
rz(2.2041722) q[1];
sx q[1];
rz(-1.2317069) q[1];
sx q[1];
rz(2.6984528) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4454058) q[0];
sx q[0];
rz(-1.6517795) q[0];
sx q[0];
rz(3.0565673) q[0];
rz(2.8117026) q[2];
sx q[2];
rz(-0.99890814) q[2];
sx q[2];
rz(-3.013333) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7305371) q[1];
sx q[1];
rz(-1.2403508) q[1];
sx q[1];
rz(-0.6544658) q[1];
rz(-pi) q[2];
rz(-1.2594873) q[3];
sx q[3];
rz(-1.7260625) q[3];
sx q[3];
rz(-0.68607611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0605269) q[2];
sx q[2];
rz(-2.492283) q[2];
sx q[2];
rz(3.013986) q[2];
rz(-0.29868948) q[3];
sx q[3];
rz(-1.9689711) q[3];
sx q[3];
rz(2.7711788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5304607) q[0];
sx q[0];
rz(-1.8831384) q[0];
sx q[0];
rz(0.25181121) q[0];
rz(0.61027377) q[1];
sx q[1];
rz(-2.2563969) q[1];
sx q[1];
rz(1.0856249) q[1];
rz(3.0771607) q[2];
sx q[2];
rz(-0.85307912) q[2];
sx q[2];
rz(-1.0882594) q[2];
rz(-1.1202624) q[3];
sx q[3];
rz(-1.7623175) q[3];
sx q[3];
rz(0.13468066) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
