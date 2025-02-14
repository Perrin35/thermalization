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
rz(-2.1433266) q[0];
sx q[0];
rz(2.7531667) q[0];
rz(1.9198963) q[1];
sx q[1];
rz(-0.95284) q[1];
sx q[1];
rz(3.067692) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.184351) q[0];
sx q[0];
rz(-1.3660407) q[0];
sx q[0];
rz(-0.53157579) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1178463) q[2];
sx q[2];
rz(-2.2226035) q[2];
sx q[2];
rz(-1.786916) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5351341) q[1];
sx q[1];
rz(-0.89790895) q[1];
sx q[1];
rz(0.4663286) q[1];
rz(-0.8624997) q[3];
sx q[3];
rz(-2.8180606) q[3];
sx q[3];
rz(2.0379553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.6834324) q[2];
sx q[2];
rz(-1.4417803) q[2];
sx q[2];
rz(-1.8140351) q[2];
rz(-0.96539998) q[3];
sx q[3];
rz(-2.5520971) q[3];
sx q[3];
rz(-2.086967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4210159) q[0];
sx q[0];
rz(-0.0076616658) q[0];
sx q[0];
rz(-1.7505919) q[0];
rz(-1.1945456) q[1];
sx q[1];
rz(-0.60940131) q[1];
sx q[1];
rz(-2.4360099) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6906889) q[0];
sx q[0];
rz(-2.2165059) q[0];
sx q[0];
rz(2.5049107) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1384355) q[2];
sx q[2];
rz(-0.92022824) q[2];
sx q[2];
rz(2.7108266) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7930562) q[1];
sx q[1];
rz(-0.67386857) q[1];
sx q[1];
rz(2.1516031) q[1];
x q[2];
rz(-0.59986214) q[3];
sx q[3];
rz(-1.9981435) q[3];
sx q[3];
rz(-0.35926705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.07448639) q[2];
sx q[2];
rz(-1.4718082) q[2];
sx q[2];
rz(2.2763695) q[2];
rz(-1.5557965) q[3];
sx q[3];
rz(-0.14388789) q[3];
sx q[3];
rz(-1.5998862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9928352) q[0];
sx q[0];
rz(-0.8388297) q[0];
sx q[0];
rz(1.6194153) q[0];
rz(1.7332227) q[1];
sx q[1];
rz(-0.38196462) q[1];
sx q[1];
rz(-0.24040374) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6344389) q[0];
sx q[0];
rz(-0.85570645) q[0];
sx q[0];
rz(-1.8832182) q[0];
rz(-2.9536134) q[2];
sx q[2];
rz(-1.9620634) q[2];
sx q[2];
rz(0.94598929) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.93546219) q[1];
sx q[1];
rz(-1.4656126) q[1];
sx q[1];
rz(1.5291924) q[1];
x q[2];
rz(1.4355074) q[3];
sx q[3];
rz(-1.093089) q[3];
sx q[3];
rz(-2.6729134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1323041) q[2];
sx q[2];
rz(-1.4762286) q[2];
sx q[2];
rz(-1.08584) q[2];
rz(-3.0911607) q[3];
sx q[3];
rz(-1.7303053) q[3];
sx q[3];
rz(-2.1850695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10629912) q[0];
sx q[0];
rz(-1.4426008) q[0];
sx q[0];
rz(-2.371149) q[0];
rz(0.92535198) q[1];
sx q[1];
rz(-1.4848361) q[1];
sx q[1];
rz(0.83522183) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2352426) q[0];
sx q[0];
rz(-0.49717227) q[0];
sx q[0];
rz(3.1072561) q[0];
rz(1.6603819) q[2];
sx q[2];
rz(-1.3302186) q[2];
sx q[2];
rz(-0.72270715) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.28047565) q[1];
sx q[1];
rz(-2.0973699) q[1];
sx q[1];
rz(-3.1269424) q[1];
x q[2];
rz(1.4553444) q[3];
sx q[3];
rz(-1.770062) q[3];
sx q[3];
rz(-1.5448567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.21021065) q[2];
sx q[2];
rz(-1.8385734) q[2];
sx q[2];
rz(2.4118928) q[2];
rz(0.99916712) q[3];
sx q[3];
rz(-1.3068643) q[3];
sx q[3];
rz(-2.6877747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5515185) q[0];
sx q[0];
rz(-2.9158264) q[0];
sx q[0];
rz(-2.1424868) q[0];
rz(-2.0935811) q[1];
sx q[1];
rz(-2.4297355) q[1];
sx q[1];
rz(-2.5968831) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22337974) q[0];
sx q[0];
rz(-1.8210016) q[0];
sx q[0];
rz(1.7556719) q[0];
rz(-pi) q[1];
rz(-0.66584058) q[2];
sx q[2];
rz(-0.84402914) q[2];
sx q[2];
rz(-0.11451463) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5156158) q[1];
sx q[1];
rz(-0.57913172) q[1];
sx q[1];
rz(2.7939151) q[1];
x q[2];
rz(-2.7361511) q[3];
sx q[3];
rz(-1.5675987) q[3];
sx q[3];
rz(-0.965626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5576632) q[2];
sx q[2];
rz(-1.9926535) q[2];
sx q[2];
rz(2.124713) q[2];
rz(2.478638) q[3];
sx q[3];
rz(-0.69279492) q[3];
sx q[3];
rz(-1.8069161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3589631) q[0];
sx q[0];
rz(-2.6435659) q[0];
sx q[0];
rz(1.9899415) q[0];
rz(-2.2478814) q[1];
sx q[1];
rz(-1.7620554) q[1];
sx q[1];
rz(-3.025257) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6457845) q[0];
sx q[0];
rz(-2.6538355) q[0];
sx q[0];
rz(1.8013823) q[0];
rz(3.1239188) q[2];
sx q[2];
rz(-1.2264614) q[2];
sx q[2];
rz(1.5196821) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.73977208) q[1];
sx q[1];
rz(-1.7171613) q[1];
sx q[1];
rz(-2.1256281) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74832423) q[3];
sx q[3];
rz(-1.6295506) q[3];
sx q[3];
rz(0.74225451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5748888) q[2];
sx q[2];
rz(-2.7723007) q[2];
sx q[2];
rz(2.134038) q[2];
rz(-1.9762074) q[3];
sx q[3];
rz(-0.27519614) q[3];
sx q[3];
rz(-0.63024855) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4646869) q[0];
sx q[0];
rz(-2.6755264) q[0];
sx q[0];
rz(-0.16978547) q[0];
rz(-2.4618497) q[1];
sx q[1];
rz(-0.83234537) q[1];
sx q[1];
rz(2.2198417) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0158495) q[0];
sx q[0];
rz(-1.4471869) q[0];
sx q[0];
rz(-2.8829734) q[0];
rz(-pi) q[1];
x q[1];
rz(2.970457) q[2];
sx q[2];
rz(-2.830392) q[2];
sx q[2];
rz(-0.61146525) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5853928) q[1];
sx q[1];
rz(-1.0784282) q[1];
sx q[1];
rz(0.3979759) q[1];
rz(-pi) q[2];
rz(2.8677651) q[3];
sx q[3];
rz(-1.3845155) q[3];
sx q[3];
rz(2.8485988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13806954) q[2];
sx q[2];
rz(-1.6602844) q[2];
sx q[2];
rz(-1.3235693) q[2];
rz(-2.0924163) q[3];
sx q[3];
rz(-1.1314355) q[3];
sx q[3];
rz(1.7607035) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84999371) q[0];
sx q[0];
rz(-2.0065362) q[0];
sx q[0];
rz(-2.3681613) q[0];
rz(-2.6749581) q[1];
sx q[1];
rz(-2.0687053) q[1];
sx q[1];
rz(-3.1007865) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.651603) q[0];
sx q[0];
rz(-1.7035653) q[0];
sx q[0];
rz(3.1285355) q[0];
rz(0.86589028) q[2];
sx q[2];
rz(-2.299832) q[2];
sx q[2];
rz(-1.2129606) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.99090138) q[1];
sx q[1];
rz(-2.5579431) q[1];
sx q[1];
rz(1.4787514) q[1];
x q[2];
rz(-1.8383305) q[3];
sx q[3];
rz(-2.5604381) q[3];
sx q[3];
rz(-2.0460956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1524973) q[2];
sx q[2];
rz(-2.5941807) q[2];
sx q[2];
rz(0.034817783) q[2];
rz(-3.1133437) q[3];
sx q[3];
rz(-1.0476799) q[3];
sx q[3];
rz(-0.69407535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5038274) q[0];
sx q[0];
rz(-0.70264188) q[0];
sx q[0];
rz(-2.0565865) q[0];
rz(-1.7833692) q[1];
sx q[1];
rz(-0.63301507) q[1];
sx q[1];
rz(1.0795275) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5732291) q[0];
sx q[0];
rz(-2.1255593) q[0];
sx q[0];
rz(2.8994096) q[0];
rz(-pi) q[1];
rz(-2.8582879) q[2];
sx q[2];
rz(-0.46479169) q[2];
sx q[2];
rz(1.9329485) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.92380556) q[1];
sx q[1];
rz(-1.3239685) q[1];
sx q[1];
rz(1.2533761) q[1];
x q[2];
rz(2.1302159) q[3];
sx q[3];
rz(-2.3314948) q[3];
sx q[3];
rz(-1.427785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8980155) q[2];
sx q[2];
rz(-1.54553) q[2];
sx q[2];
rz(-3.0698981) q[2];
rz(-0.60595766) q[3];
sx q[3];
rz(-2.3809483) q[3];
sx q[3];
rz(3.0683556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3074985) q[0];
sx q[0];
rz(-1.5102757) q[0];
sx q[0];
rz(0.50258762) q[0];
rz(1.7723627) q[1];
sx q[1];
rz(-1.2255729) q[1];
sx q[1];
rz(-2.9659081) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0640819) q[0];
sx q[0];
rz(-0.33987576) q[0];
sx q[0];
rz(1.3114291) q[0];
rz(2.2403342) q[2];
sx q[2];
rz(-0.84543919) q[2];
sx q[2];
rz(-2.8536316) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8524234) q[1];
sx q[1];
rz(-0.90327493) q[1];
sx q[1];
rz(-0.65594419) q[1];
rz(-pi) q[2];
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
x q[1];
rz(2.6548369) q[2];
sx q[2];
rz(-1.2590057) q[2];
sx q[2];
rz(-0.34461018) q[2];
rz(1.9136072) q[3];
sx q[3];
rz(-0.69209185) q[3];
sx q[3];
rz(-0.96755782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6268613) q[0];
sx q[0];
rz(-1.3118962) q[0];
sx q[0];
rz(0.94919039) q[0];
rz(1.3728036) q[1];
sx q[1];
rz(-2.1887442) q[1];
sx q[1];
rz(-1.8123117) q[1];
rz(-2.4288154) q[2];
sx q[2];
rz(-1.3679124) q[2];
sx q[2];
rz(-1.5001679) q[2];
rz(2.9459841) q[3];
sx q[3];
rz(-0.8630639) q[3];
sx q[3];
rz(0.56474781) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
