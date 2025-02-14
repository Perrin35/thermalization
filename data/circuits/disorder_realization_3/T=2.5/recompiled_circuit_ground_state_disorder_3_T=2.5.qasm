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
rz(-1.87529110908508) q[0];
sx q[0];
rz(3.54850617249543) q[0];
sx q[0];
rz(9.22235280870601) q[0];
rz(1.83167326450348) q[1];
sx q[1];
rz(4.28694215615327) q[1];
sx q[1];
rz(11.5776791334073) q[1];
cx q[1],q[0];
rz(1.76449728012085) q[0];
sx q[0];
rz(5.45410791237886) q[0];
sx q[0];
rz(8.11202869414493) q[0];
rz(2.62767720222473) q[2];
sx q[2];
rz(2.98364532192285) q[2];
sx q[2];
rz(5.22280166148349) q[2];
cx q[2],q[1];
rz(-0.0844170153141022) q[1];
sx q[1];
rz(4.93817725976045) q[1];
sx q[1];
rz(8.79762235879108) q[1];
rz(-0.268817484378815) q[3];
sx q[3];
rz(3.58146661718423) q[3];
sx q[3];
rz(8.30379590987369) q[3];
cx q[3],q[2];
rz(0.536283552646637) q[2];
sx q[2];
rz(3.94318959315354) q[2];
sx q[2];
rz(10.4762001991193) q[2];
rz(-1.50315225124359) q[3];
sx q[3];
rz(4.5548541863733) q[3];
sx q[3];
rz(9.66004156171485) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.19920921325684) q[0];
sx q[0];
rz(2.11140445073182) q[0];
sx q[0];
rz(9.69879404305621) q[0];
rz(0.383852988481522) q[1];
sx q[1];
rz(3.77420649130876) q[1];
sx q[1];
rz(10.7455597877423) q[1];
cx q[1],q[0];
rz(0.216525539755821) q[0];
sx q[0];
rz(4.55379739602143) q[0];
sx q[0];
rz(10.9301699161451) q[0];
rz(-1.97813773155212) q[2];
sx q[2];
rz(3.67766723235185) q[2];
sx q[2];
rz(11.6323253869931) q[2];
cx q[2],q[1];
rz(2.74104237556458) q[1];
sx q[1];
rz(2.29307511647279) q[1];
sx q[1];
rz(10.9321731090467) q[1];
rz(-0.781353235244751) q[3];
sx q[3];
rz(2.20092103083665) q[3];
sx q[3];
rz(10.3321561574857) q[3];
cx q[3],q[2];
rz(2.54478335380554) q[2];
sx q[2];
rz(4.88045981724794) q[2];
sx q[2];
rz(10.8242212295453) q[2];
rz(0.607473731040955) q[3];
sx q[3];
rz(1.67392686207826) q[3];
sx q[3];
rz(9.81529954671069) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.529432773590088) q[0];
sx q[0];
rz(2.2859891970926) q[0];
sx q[0];
rz(9.10554719566509) q[0];
rz(0.0514797382056713) q[1];
sx q[1];
rz(5.10365703900392) q[1];
sx q[1];
rz(7.33974406718417) q[1];
cx q[1],q[0];
rz(-0.385045439004898) q[0];
sx q[0];
rz(2.96998848219449) q[0];
sx q[0];
rz(9.11848042010471) q[0];
rz(-0.705510675907135) q[2];
sx q[2];
rz(4.74924627144868) q[2];
sx q[2];
rz(11.2591801643293) q[2];
cx q[2],q[1];
rz(-0.332568615674973) q[1];
sx q[1];
rz(4.5377301295572) q[1];
sx q[1];
rz(10.0578199982564) q[1];
rz(1.66328358650208) q[3];
sx q[3];
rz(4.31068721612031) q[3];
sx q[3];
rz(10.1432794690053) q[3];
cx q[3],q[2];
rz(-1.94873785972595) q[2];
sx q[2];
rz(4.57338348229463) q[2];
sx q[2];
rz(10.7693271398465) q[2];
rz(0.917908489704132) q[3];
sx q[3];
rz(2.50360056956346) q[3];
sx q[3];
rz(11.3105451822202) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.464808791875839) q[0];
sx q[0];
rz(3.73216125567491) q[0];
sx q[0];
rz(11.0574465751569) q[0];
rz(-0.50679874420166) q[1];
sx q[1];
rz(4.3549965937906) q[1];
sx q[1];
rz(9.81533617376491) q[1];
cx q[1],q[0];
rz(1.20295011997223) q[0];
sx q[0];
rz(3.14999860537285) q[0];
sx q[0];
rz(10.4672659397046) q[0];
rz(0.190436810255051) q[2];
sx q[2];
rz(1.94452968438203) q[2];
sx q[2];
rz(8.44294432400867) q[2];
cx q[2],q[1];
rz(1.64667546749115) q[1];
sx q[1];
rz(3.7611764391237) q[1];
sx q[1];
rz(9.47365367635294) q[1];
rz(2.08708763122559) q[3];
sx q[3];
rz(3.32328739960725) q[3];
sx q[3];
rz(8.80392280816242) q[3];
cx q[3],q[2];
rz(-1.89166963100433) q[2];
sx q[2];
rz(2.21077832778031) q[2];
sx q[2];
rz(11.6699731111447) q[2];
rz(-0.916263818740845) q[3];
sx q[3];
rz(4.07742676337297) q[3];
sx q[3];
rz(10.385522401325) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.591782808303833) q[0];
sx q[0];
rz(3.49742484291131) q[0];
sx q[0];
rz(9.89497462510272) q[0];
rz(-1.40379405021667) q[1];
sx q[1];
rz(3.84198847611482) q[1];
sx q[1];
rz(7.51087770461246) q[1];
cx q[1],q[0];
rz(-1.02253425121307) q[0];
sx q[0];
rz(4.7187489589029) q[0];
sx q[0];
rz(12.2962591409604) q[0];
rz(2.71409845352173) q[2];
sx q[2];
rz(4.20133963425691) q[2];
sx q[2];
rz(11.0781678914945) q[2];
cx q[2],q[1];
rz(-1.22540819644928) q[1];
sx q[1];
rz(5.63403669198091) q[1];
sx q[1];
rz(9.29937360285922) q[1];
rz(0.310570299625397) q[3];
sx q[3];
rz(5.22938457329805) q[3];
sx q[3];
rz(9.65111940204307) q[3];
cx q[3],q[2];
rz(0.970518350601196) q[2];
sx q[2];
rz(3.85434887011582) q[2];
sx q[2];
rz(10.8090324163358) q[2];
rz(0.61128830909729) q[3];
sx q[3];
rz(4.52117410500581) q[3];
sx q[3];
rz(9.02864170669719) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.14320814609528) q[0];
sx q[0];
rz(5.27135268052156) q[0];
sx q[0];
rz(9.97692528962299) q[0];
rz(2.41190075874329) q[1];
sx q[1];
rz(2.15398600895936) q[1];
sx q[1];
rz(8.42302033900424) q[1];
cx q[1],q[0];
rz(-0.692437887191772) q[0];
sx q[0];
rz(4.57913378079469) q[0];
sx q[0];
rz(7.93521616458102) q[0];
rz(-0.628436267375946) q[2];
sx q[2];
rz(5.53714242775971) q[2];
sx q[2];
rz(7.50682685374423) q[2];
cx q[2],q[1];
rz(-0.616484701633453) q[1];
sx q[1];
rz(4.89614728291566) q[1];
sx q[1];
rz(10.2241185069005) q[1];
rz(0.891734421253204) q[3];
sx q[3];
rz(1.15213945706422) q[3];
sx q[3];
rz(8.4763747215192) q[3];
cx q[3],q[2];
rz(0.457276493310928) q[2];
sx q[2];
rz(4.49068990548188) q[2];
sx q[2];
rz(11.2229972839276) q[2];
rz(0.724304795265198) q[3];
sx q[3];
rz(1.93800059159333) q[3];
sx q[3];
rz(9.12211380004092) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.901720762252808) q[0];
sx q[0];
rz(4.47387495835359) q[0];
sx q[0];
rz(10.1501943826596) q[0];
rz(-1.8431009054184) q[1];
sx q[1];
rz(3.59386819799478) q[1];
sx q[1];
rz(12.287669157974) q[1];
cx q[1],q[0];
rz(0.505100786685944) q[0];
sx q[0];
rz(4.17676285107667) q[0];
sx q[0];
rz(9.24779916404887) q[0];
rz(-2.45298552513123) q[2];
sx q[2];
rz(5.7857028563791) q[2];
sx q[2];
rz(8.90537527798816) q[2];
cx q[2],q[1];
rz(0.781657516956329) q[1];
sx q[1];
rz(5.0167374928766) q[1];
sx q[1];
rz(11.8370022535245) q[1];
rz(-0.549460649490356) q[3];
sx q[3];
rz(3.96449604828889) q[3];
sx q[3];
rz(11.3149814367215) q[3];
cx q[3],q[2];
rz(0.384103864431381) q[2];
sx q[2];
rz(3.40833065112168) q[2];
sx q[2];
rz(8.43662992715045) q[2];
rz(3.24654769897461) q[3];
sx q[3];
rz(4.98493650754029) q[3];
sx q[3];
rz(9.65611434578105) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.688954830169678) q[0];
sx q[0];
rz(5.48204842408235) q[0];
sx q[0];
rz(8.80538586377307) q[0];
rz(3.15648484230042) q[1];
sx q[1];
rz(4.03193286259706) q[1];
sx q[1];
rz(6.96488616465732) q[1];
cx q[1],q[0];
rz(1.83121633529663) q[0];
sx q[0];
rz(4.99995294411714) q[0];
sx q[0];
rz(9.98470768927737) q[0];
rz(-3.24370455741882) q[2];
sx q[2];
rz(3.80925491650636) q[2];
sx q[2];
rz(13.4073996305387) q[2];
cx q[2],q[1];
rz(1.35779631137848) q[1];
sx q[1];
rz(5.06397691567475) q[1];
sx q[1];
rz(9.43200694014459) q[1];
rz(0.210582450032234) q[3];
sx q[3];
rz(3.89040050108964) q[3];
sx q[3];
rz(10.5960877895276) q[3];
cx q[3],q[2];
rz(0.909984886646271) q[2];
sx q[2];
rz(3.94474253256852) q[2];
sx q[2];
rz(12.226713871948) q[2];
rz(0.613994419574738) q[3];
sx q[3];
rz(4.52117052872712) q[3];
sx q[3];
rz(10.9865034580152) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.15641076862812) q[0];
sx q[0];
rz(3.48973202903802) q[0];
sx q[0];
rz(8.57038477658435) q[0];
rz(0.593354523181915) q[1];
sx q[1];
rz(4.17363849480683) q[1];
sx q[1];
rz(9.72012231349155) q[1];
cx q[1],q[0];
rz(-0.127960622310638) q[0];
sx q[0];
rz(4.67951503594453) q[0];
sx q[0];
rz(11.0508782625119) q[0];
rz(-0.0507120564579964) q[2];
sx q[2];
rz(3.5598220547014) q[2];
sx q[2];
rz(9.02755332588359) q[2];
cx q[2],q[1];
rz(2.01573896408081) q[1];
sx q[1];
rz(5.39120403130586) q[1];
sx q[1];
rz(8.48292622565433) q[1];
rz(2.80286502838135) q[3];
sx q[3];
rz(2.81451878150041) q[3];
sx q[3];
rz(7.33363745211765) q[3];
cx q[3],q[2];
rz(-2.97965431213379) q[2];
sx q[2];
rz(4.82478943665559) q[2];
sx q[2];
rz(11.8845386266629) q[2];
rz(-1.94090592861176) q[3];
sx q[3];
rz(3.92939076026017) q[3];
sx q[3];
rz(9.86062059401675) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.110726304352283) q[0];
sx q[0];
rz(4.42152813275392) q[0];
sx q[0];
rz(9.93064138888522) q[0];
rz(1.01816308498383) q[1];
sx q[1];
rz(4.84804263909394) q[1];
sx q[1];
rz(8.56431219576999) q[1];
cx q[1],q[0];
rz(1.84485256671906) q[0];
sx q[0];
rz(0.895497949915477) q[0];
sx q[0];
rz(9.3642231658022) q[0];
rz(0.292975962162018) q[2];
sx q[2];
rz(4.57618406613404) q[2];
sx q[2];
rz(9.95538440941974) q[2];
cx q[2],q[1];
rz(-0.0827493593096733) q[1];
sx q[1];
rz(2.89020210702951) q[1];
sx q[1];
rz(10.8631649970929) q[1];
rz(1.00603723526001) q[3];
sx q[3];
rz(4.1755132993036) q[3];
sx q[3];
rz(10.3443245053212) q[3];
cx q[3],q[2];
rz(1.08184170722961) q[2];
sx q[2];
rz(1.95307186444337) q[2];
sx q[2];
rz(8.84158310889407) q[2];
rz(0.0010631337063387) q[3];
sx q[3];
rz(4.07261529763276) q[3];
sx q[3];
rz(9.91812700628444) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.902768313884735) q[0];
sx q[0];
rz(2.45841631491716) q[0];
sx q[0];
rz(9.66658035515949) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(2.52316641807556) q[1];
sx q[1];
rz(2.89635878999765) q[1];
sx q[1];
rz(9.07811582683727) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(0.990509033203125) q[2];
sx q[2];
rz(3.60975903471047) q[2];
sx q[2];
rz(12.0512196779172) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(1.74771523475647) q[3];
sx q[3];
rz(4.78907218773896) q[3];
sx q[3];
rz(8.44489959477588) q[3];
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
