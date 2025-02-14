OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5724343) q[0];
sx q[0];
rz(7.2814514) q[0];
sx q[0];
rz(9.8132039) q[0];
rz(-1.2216964) q[1];
sx q[1];
rz(-2.1887527) q[1];
sx q[1];
rz(0.0739007) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6361481) q[0];
sx q[0];
rz(-1.0514604) q[0];
sx q[0];
rz(1.8072007) q[0];
rz(-pi) q[1];
rz(1.6019033) q[2];
sx q[2];
rz(-2.4894161) q[2];
sx q[2];
rz(1.3155441) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5351341) q[1];
sx q[1];
rz(-0.89790895) q[1];
sx q[1];
rz(0.4663286) q[1];
x q[2];
rz(0.8624997) q[3];
sx q[3];
rz(-0.32353208) q[3];
sx q[3];
rz(-1.1036373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4581603) q[2];
sx q[2];
rz(-1.4417803) q[2];
sx q[2];
rz(1.8140351) q[2];
rz(-0.96539998) q[3];
sx q[3];
rz(-2.5520971) q[3];
sx q[3];
rz(-2.086967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7205768) q[0];
sx q[0];
rz(-0.0076616658) q[0];
sx q[0];
rz(1.3910008) q[0];
rz(1.1945456) q[1];
sx q[1];
rz(-0.60940131) q[1];
sx q[1];
rz(2.4360099) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53853453) q[0];
sx q[0];
rz(-1.0760191) q[0];
sx q[0];
rz(2.323708) q[0];
rz(2.6383361) q[2];
sx q[2];
rz(-2.3782298) q[2];
sx q[2];
rz(-2.9212375) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3485364) q[1];
sx q[1];
rz(-0.67386857) q[1];
sx q[1];
rz(-0.98998951) q[1];
x q[2];
rz(0.67882244) q[3];
sx q[3];
rz(-2.4206483) q[3];
sx q[3];
rz(-1.7563422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0671063) q[2];
sx q[2];
rz(-1.6697845) q[2];
sx q[2];
rz(2.2763695) q[2];
rz(-1.5857961) q[3];
sx q[3];
rz(-0.14388789) q[3];
sx q[3];
rz(-1.5417064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9928352) q[0];
sx q[0];
rz(-2.302763) q[0];
sx q[0];
rz(1.5221773) q[0];
rz(-1.4083699) q[1];
sx q[1];
rz(-2.759628) q[1];
sx q[1];
rz(0.24040374) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8692476) q[0];
sx q[0];
rz(-1.8050021) q[0];
sx q[0];
rz(-0.73973013) q[0];
rz(-pi) q[1];
rz(1.1454606) q[2];
sx q[2];
rz(-0.43197235) q[2];
sx q[2];
rz(-1.4086823) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2061305) q[1];
sx q[1];
rz(-1.67598) q[1];
sx q[1];
rz(1.6124003) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4355074) q[3];
sx q[3];
rz(-1.093089) q[3];
sx q[3];
rz(2.6729134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1323041) q[2];
sx q[2];
rz(-1.4762286) q[2];
sx q[2];
rz(-1.08584) q[2];
rz(-0.050431937) q[3];
sx q[3];
rz(-1.4112873) q[3];
sx q[3];
rz(0.95652318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0352935) q[0];
sx q[0];
rz(-1.6989919) q[0];
sx q[0];
rz(-0.77044368) q[0];
rz(-0.92535198) q[1];
sx q[1];
rz(-1.4848361) q[1];
sx q[1];
rz(2.3063708) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90635008) q[0];
sx q[0];
rz(-0.49717227) q[0];
sx q[0];
rz(-0.034336523) q[0];
rz(-0.24150924) q[2];
sx q[2];
rz(-1.4837974) q[2];
sx q[2];
rz(-2.3149025) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2829574) q[1];
sx q[1];
rz(-1.5581308) q[1];
sx q[1];
rz(1.0441761) q[1];
x q[2];
rz(-1.4553444) q[3];
sx q[3];
rz(-1.3715306) q[3];
sx q[3];
rz(-1.5448567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.21021065) q[2];
sx q[2];
rz(-1.3030193) q[2];
sx q[2];
rz(-2.4118928) q[2];
rz(0.99916712) q[3];
sx q[3];
rz(-1.3068643) q[3];
sx q[3];
rz(-2.6877747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5900742) q[0];
sx q[0];
rz(-2.9158264) q[0];
sx q[0];
rz(0.99910587) q[0];
rz(2.0935811) q[1];
sx q[1];
rz(-0.71185714) q[1];
sx q[1];
rz(-2.5968831) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9182129) q[0];
sx q[0];
rz(-1.8210016) q[0];
sx q[0];
rz(1.7556719) q[0];
rz(-2.1779973) q[2];
sx q[2];
rz(-0.94258833) q[2];
sx q[2];
rz(2.387418) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0346086) q[1];
sx q[1];
rz(-1.0303127) q[1];
sx q[1];
rz(1.3515737) q[1];
rz(-pi) q[2];
rz(0.0081069907) q[3];
sx q[3];
rz(-2.7361392) q[3];
sx q[3];
rz(2.5289726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5576632) q[2];
sx q[2];
rz(-1.9926535) q[2];
sx q[2];
rz(-2.124713) q[2];
rz(0.66295463) q[3];
sx q[3];
rz(-0.69279492) q[3];
sx q[3];
rz(1.8069161) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7826295) q[0];
sx q[0];
rz(-0.49802676) q[0];
sx q[0];
rz(1.9899415) q[0];
rz(0.89371124) q[1];
sx q[1];
rz(-1.3795373) q[1];
sx q[1];
rz(3.025257) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75555629) q[0];
sx q[0];
rz(-2.0445637) q[0];
sx q[0];
rz(-3.0209345) q[0];
x q[1];
rz(-1.5215559) q[2];
sx q[2];
rz(-0.34477012) q[2];
sx q[2];
rz(1.5719959) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.92116881) q[1];
sx q[1];
rz(-1.0225778) q[1];
sx q[1];
rz(2.9698644) q[1];
x q[2];
rz(-0.086236091) q[3];
sx q[3];
rz(-2.3914118) q[3];
sx q[3];
rz(2.3762109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5748888) q[2];
sx q[2];
rz(-0.36929193) q[2];
sx q[2];
rz(-2.134038) q[2];
rz(1.1653853) q[3];
sx q[3];
rz(-0.27519614) q[3];
sx q[3];
rz(-0.63024855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6769058) q[0];
sx q[0];
rz(-2.6755264) q[0];
sx q[0];
rz(0.16978547) q[0];
rz(-2.4618497) q[1];
sx q[1];
rz(-2.3092473) q[1];
sx q[1];
rz(-2.2198417) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1257432) q[0];
sx q[0];
rz(-1.4471869) q[0];
sx q[0];
rz(-2.8829734) q[0];
x q[1];
rz(-2.970457) q[2];
sx q[2];
rz(-2.830392) q[2];
sx q[2];
rz(0.61146525) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8584916) q[1];
sx q[1];
rz(-2.518961) q[1];
sx q[1];
rz(2.1964392) q[1];
x q[2];
rz(0.60867057) q[3];
sx q[3];
rz(-2.8117315) q[3];
sx q[3];
rz(0.69486991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0035231) q[2];
sx q[2];
rz(-1.4813083) q[2];
sx q[2];
rz(1.3235693) q[2];
rz(-2.0924163) q[3];
sx q[3];
rz(-2.0101571) q[3];
sx q[3];
rz(1.3808892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2915989) q[0];
sx q[0];
rz(-1.1350564) q[0];
sx q[0];
rz(2.3681613) q[0];
rz(-0.4666346) q[1];
sx q[1];
rz(-2.0687053) q[1];
sx q[1];
rz(-0.040806142) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.651603) q[0];
sx q[0];
rz(-1.7035653) q[0];
sx q[0];
rz(-0.013057166) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62760533) q[2];
sx q[2];
rz(-0.96663953) q[2];
sx q[2];
rz(-1.0221572) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0405214) q[1];
sx q[1];
rz(-2.1516529) q[1];
sx q[1];
rz(3.0809666) q[1];
x q[2];
rz(-1.3032622) q[3];
sx q[3];
rz(-2.5604381) q[3];
sx q[3];
rz(-1.095497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9890954) q[2];
sx q[2];
rz(-2.5941807) q[2];
sx q[2];
rz(-3.1067749) q[2];
rz(3.1133437) q[3];
sx q[3];
rz(-1.0476799) q[3];
sx q[3];
rz(-2.4475173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5038274) q[0];
sx q[0];
rz(-2.4389508) q[0];
sx q[0];
rz(-1.0850061) q[0];
rz(-1.7833692) q[1];
sx q[1];
rz(-2.5085776) q[1];
sx q[1];
rz(2.0620652) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5732291) q[0];
sx q[0];
rz(-1.0160334) q[0];
sx q[0];
rz(-2.8994096) q[0];
rz(-pi) q[1];
rz(-0.28330477) q[2];
sx q[2];
rz(-0.46479169) q[2];
sx q[2];
rz(-1.9329485) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.72708817) q[1];
sx q[1];
rz(-1.263317) q[1];
sx q[1];
rz(-2.8823463) q[1];
x q[2];
rz(-2.2983389) q[3];
sx q[3];
rz(-1.9653659) q[3];
sx q[3];
rz(-0.55052557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8980155) q[2];
sx q[2];
rz(-1.5960627) q[2];
sx q[2];
rz(-0.071694516) q[2];
rz(0.60595766) q[3];
sx q[3];
rz(-0.76064435) q[3];
sx q[3];
rz(-0.073237091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3074985) q[0];
sx q[0];
rz(-1.5102757) q[0];
sx q[0];
rz(2.639005) q[0];
rz(-1.7723627) q[1];
sx q[1];
rz(-1.9160198) q[1];
sx q[1];
rz(0.1756846) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7384199) q[0];
sx q[0];
rz(-1.4851928) q[0];
sx q[0];
rz(-1.241472) q[0];
rz(-0.90125842) q[2];
sx q[2];
rz(-2.2961535) q[2];
sx q[2];
rz(-0.28796107) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.726242) q[1];
sx q[1];
rz(-1.0712855) q[1];
sx q[1];
rz(-2.3535092) q[1];
x q[2];
rz(-2.7501372) q[3];
sx q[3];
rz(-2.1717697) q[3];
sx q[3];
rz(1.6319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6548369) q[2];
sx q[2];
rz(-1.882587) q[2];
sx q[2];
rz(2.7969825) q[2];
rz(-1.2279855) q[3];
sx q[3];
rz(-2.4495008) q[3];
sx q[3];
rz(-2.1740348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5147314) q[0];
sx q[0];
rz(-1.3118962) q[0];
sx q[0];
rz(0.94919039) q[0];
rz(-1.3728036) q[1];
sx q[1];
rz(-0.95284843) q[1];
sx q[1];
rz(1.3292809) q[1];
rz(-2.8368159) q[2];
sx q[2];
rz(-0.73619107) q[2];
sx q[2];
rz(0.29965055) q[2];
rz(-1.7941734) q[3];
sx q[3];
rz(-2.4118524) q[3];
sx q[3];
rz(0.86058544) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
