OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.9450077) q[0];
sx q[0];
rz(-0.44214806) q[0];
sx q[0];
rz(-2.2801939) q[0];
rz(-1.6317033) q[1];
sx q[1];
rz(-1.3246526) q[1];
sx q[1];
rz(-3.1060001) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0738348) q[0];
sx q[0];
rz(-1.375421) q[0];
sx q[0];
rz(-2.8003576) q[0];
x q[1];
rz(-2.8350949) q[2];
sx q[2];
rz(-1.674429) q[2];
sx q[2];
rz(0.95970861) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5000677) q[1];
sx q[1];
rz(-1.6065805) q[1];
sx q[1];
rz(-0.72662093) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5930327) q[3];
sx q[3];
rz(-2.4614868) q[3];
sx q[3];
rz(-0.82397991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5888136) q[2];
sx q[2];
rz(-1.3014883) q[2];
sx q[2];
rz(2.0471052) q[2];
rz(-1.6257446) q[3];
sx q[3];
rz(-2.2690014) q[3];
sx q[3];
rz(-2.998013) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0482408) q[0];
sx q[0];
rz(-2.1513262) q[0];
sx q[0];
rz(-2.7675203) q[0];
rz(-2.2834942) q[1];
sx q[1];
rz(-0.94081196) q[1];
sx q[1];
rz(-2.0657952) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6158524) q[0];
sx q[0];
rz(-1.6739556) q[0];
sx q[0];
rz(3.0391058) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38041842) q[2];
sx q[2];
rz(-2.5606692) q[2];
sx q[2];
rz(0.38512938) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3871284) q[1];
sx q[1];
rz(-0.417388) q[1];
sx q[1];
rz(-2.8512958) q[1];
rz(-2.2298563) q[3];
sx q[3];
rz(-2.1359332) q[3];
sx q[3];
rz(-2.679058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6263803) q[2];
sx q[2];
rz(-2.7084646) q[2];
sx q[2];
rz(1.0832146) q[2];
rz(-1.8488688) q[3];
sx q[3];
rz(-2.0079398) q[3];
sx q[3];
rz(-0.92866549) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.885963) q[0];
sx q[0];
rz(-2.3704381) q[0];
sx q[0];
rz(2.6166925) q[0];
rz(1.4595754) q[1];
sx q[1];
rz(-0.73760215) q[1];
sx q[1];
rz(2.9958013) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73506671) q[0];
sx q[0];
rz(-0.13938381) q[0];
sx q[0];
rz(0.71203072) q[0];
rz(-0.86051382) q[2];
sx q[2];
rz(-2.4053221) q[2];
sx q[2];
rz(-2.3585573) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5151095) q[1];
sx q[1];
rz(-1.4806387) q[1];
sx q[1];
rz(0.43403352) q[1];
rz(2.5691312) q[3];
sx q[3];
rz(-2.2160071) q[3];
sx q[3];
rz(1.087134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8968481) q[2];
sx q[2];
rz(-1.9795828) q[2];
sx q[2];
rz(0.16540089) q[2];
rz(2.2798955) q[3];
sx q[3];
rz(-0.69178897) q[3];
sx q[3];
rz(0.19680463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8621314) q[0];
sx q[0];
rz(-2.617351) q[0];
sx q[0];
rz(-2.4133546) q[0];
rz(-0.32314745) q[1];
sx q[1];
rz(-1.0341945) q[1];
sx q[1];
rz(-1.039215) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.090492608) q[0];
sx q[0];
rz(-1.9673825) q[0];
sx q[0];
rz(1.1086247) q[0];
rz(2.2065998) q[2];
sx q[2];
rz(-1.3424917) q[2];
sx q[2];
rz(1.1479957) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7936662) q[1];
sx q[1];
rz(-1.4783585) q[1];
sx q[1];
rz(2.6562009) q[1];
rz(-0.83926852) q[3];
sx q[3];
rz(-2.2429562) q[3];
sx q[3];
rz(-1.9306926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6497961) q[2];
sx q[2];
rz(-2.0911262) q[2];
sx q[2];
rz(0.51898471) q[2];
rz(1.0745878) q[3];
sx q[3];
rz(-1.855775) q[3];
sx q[3];
rz(-0.49526596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.062926) q[0];
sx q[0];
rz(-0.92200297) q[0];
sx q[0];
rz(-2.9700188) q[0];
rz(-1.0942787) q[1];
sx q[1];
rz(-1.8405874) q[1];
sx q[1];
rz(-2.9069854) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5253403) q[0];
sx q[0];
rz(-1.8484637) q[0];
sx q[0];
rz(0.13292952) q[0];
rz(-pi) q[1];
rz(1.5182839) q[2];
sx q[2];
rz(-1.0641452) q[2];
sx q[2];
rz(-1.6135297) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7947096) q[1];
sx q[1];
rz(-0.88635072) q[1];
sx q[1];
rz(-2.3331932) q[1];
rz(-0.71013807) q[3];
sx q[3];
rz(-1.8605068) q[3];
sx q[3];
rz(-0.75677815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.91308633) q[2];
sx q[2];
rz(-1.9501016) q[2];
sx q[2];
rz(1.7388434) q[2];
rz(-2.5168822) q[3];
sx q[3];
rz(-2.1731845) q[3];
sx q[3];
rz(-1.3431842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45873555) q[0];
sx q[0];
rz(-3.1338437) q[0];
sx q[0];
rz(0.32753456) q[0];
rz(3.1134743) q[1];
sx q[1];
rz(-1.997812) q[1];
sx q[1];
rz(-0.99172529) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92627996) q[0];
sx q[0];
rz(-0.87020391) q[0];
sx q[0];
rz(-2.451137) q[0];
x q[1];
rz(2.9155511) q[2];
sx q[2];
rz(-1.4648572) q[2];
sx q[2];
rz(2.6424579) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.36054143) q[1];
sx q[1];
rz(-2.326093) q[1];
sx q[1];
rz(-0.057572854) q[1];
rz(-pi) q[2];
rz(1.5981653) q[3];
sx q[3];
rz(-2.0413766) q[3];
sx q[3];
rz(-1.1740299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2345387) q[2];
sx q[2];
rz(-1.3776642) q[2];
sx q[2];
rz(-1.6592337) q[2];
rz(1.0659069) q[3];
sx q[3];
rz(-1.9612471) q[3];
sx q[3];
rz(-2.0445686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7557573) q[0];
sx q[0];
rz(-0.11180728) q[0];
sx q[0];
rz(-0.29543153) q[0];
rz(-0.23751986) q[1];
sx q[1];
rz(-0.93370456) q[1];
sx q[1];
rz(-2.3942153) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6866561) q[0];
sx q[0];
rz(-1.0191518) q[0];
sx q[0];
rz(-2.174189) q[0];
x q[1];
rz(0.68589034) q[2];
sx q[2];
rz(-2.389156) q[2];
sx q[2];
rz(0.28960944) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.805757) q[1];
sx q[1];
rz(-2.0621967) q[1];
sx q[1];
rz(-2.7193428) q[1];
rz(-1.9456057) q[3];
sx q[3];
rz(-2.0829933) q[3];
sx q[3];
rz(0.95694816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6682917) q[2];
sx q[2];
rz(-1.1857727) q[2];
sx q[2];
rz(-2.8866923) q[2];
rz(0.21236803) q[3];
sx q[3];
rz(-2.2771211) q[3];
sx q[3];
rz(0.00096850639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2905228) q[0];
sx q[0];
rz(-1.2238598) q[0];
sx q[0];
rz(3.0176924) q[0];
rz(-0.87521416) q[1];
sx q[1];
rz(-2.3085935) q[1];
sx q[1];
rz(2.5828054) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0103115) q[0];
sx q[0];
rz(-1.049066) q[0];
sx q[0];
rz(2.8521796) q[0];
rz(-0.28382878) q[2];
sx q[2];
rz(-2.0528194) q[2];
sx q[2];
rz(1.986077) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9001101) q[1];
sx q[1];
rz(-2.0962976) q[1];
sx q[1];
rz(1.2578169) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93620091) q[3];
sx q[3];
rz(-1.9125328) q[3];
sx q[3];
rz(0.5638916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0957886) q[2];
sx q[2];
rz(-2.6498821) q[2];
sx q[2];
rz(2.4208505) q[2];
rz(0.36302429) q[3];
sx q[3];
rz(-1.2237153) q[3];
sx q[3];
rz(1.5117517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95661288) q[0];
sx q[0];
rz(-0.19421254) q[0];
sx q[0];
rz(-3.0525364) q[0];
rz(2.7748499) q[1];
sx q[1];
rz(-0.9042424) q[1];
sx q[1];
rz(-1.4595703) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90208611) q[0];
sx q[0];
rz(-1.7988822) q[0];
sx q[0];
rz(0.93241623) q[0];
rz(2.1543571) q[2];
sx q[2];
rz(-2.1676807) q[2];
sx q[2];
rz(1.0839628) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.0017346026) q[1];
sx q[1];
rz(-2.325104) q[1];
sx q[1];
rz(-1.0553318) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59883786) q[3];
sx q[3];
rz(-1.9396848) q[3];
sx q[3];
rz(1.2303331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9350932) q[2];
sx q[2];
rz(-1.3672028) q[2];
sx q[2];
rz(1.8915141) q[2];
rz(-1.9153204) q[3];
sx q[3];
rz(-2.5076702) q[3];
sx q[3];
rz(0.63898501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6298237) q[0];
sx q[0];
rz(-1.1331695) q[0];
sx q[0];
rz(1.1400219) q[0];
rz(0.58311588) q[1];
sx q[1];
rz(-1.4864328) q[1];
sx q[1];
rz(0.66782943) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7191737) q[0];
sx q[0];
rz(-1.8961217) q[0];
sx q[0];
rz(-2.3803985) q[0];
x q[1];
rz(-0.23867757) q[2];
sx q[2];
rz(-2.2584923) q[2];
sx q[2];
rz(2.4872045) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.17396233) q[1];
sx q[1];
rz(-2.6625427) q[1];
sx q[1];
rz(1.1381989) q[1];
rz(-pi) q[2];
x q[2];
rz(0.45248453) q[3];
sx q[3];
rz(-2.2960032) q[3];
sx q[3];
rz(0.11115995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2976611) q[2];
sx q[2];
rz(-0.49489489) q[2];
sx q[2];
rz(-2.3893791) q[2];
rz(0.080605896) q[3];
sx q[3];
rz(-1.3496496) q[3];
sx q[3];
rz(0.54752553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3362296) q[0];
sx q[0];
rz(-1.6148051) q[0];
sx q[0];
rz(-0.057407277) q[0];
rz(-0.84857955) q[1];
sx q[1];
rz(-2.4129557) q[1];
sx q[1];
rz(1.6853263) q[1];
rz(-2.2304429) q[2];
sx q[2];
rz(-1.0187889) q[2];
sx q[2];
rz(2.7762882) q[2];
rz(0.6497339) q[3];
sx q[3];
rz(-2.5841373) q[3];
sx q[3];
rz(2.0852603) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
