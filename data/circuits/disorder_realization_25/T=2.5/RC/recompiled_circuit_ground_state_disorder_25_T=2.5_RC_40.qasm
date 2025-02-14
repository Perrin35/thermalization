OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9659757) q[0];
sx q[0];
rz(-1.5953925) q[0];
sx q[0];
rz(1.6311128) q[0];
rz(1.0881967) q[1];
sx q[1];
rz(-1.3047941) q[1];
sx q[1];
rz(2.3545177) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2505456) q[0];
sx q[0];
rz(-2.4536687) q[0];
sx q[0];
rz(-0.1767747) q[0];
rz(-1.5163312) q[2];
sx q[2];
rz(-1.174841) q[2];
sx q[2];
rz(1.2690282) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0144498) q[1];
sx q[1];
rz(-2.9269955) q[1];
sx q[1];
rz(0.17226528) q[1];
rz(0.98104279) q[3];
sx q[3];
rz(-1.0372122) q[3];
sx q[3];
rz(2.5887845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2934908) q[2];
sx q[2];
rz(-0.11152554) q[2];
sx q[2];
rz(-0.21609406) q[2];
rz(0.50348336) q[3];
sx q[3];
rz(-1.2516021) q[3];
sx q[3];
rz(-3.0751198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1450495) q[0];
sx q[0];
rz(-1.1798877) q[0];
sx q[0];
rz(-2.766372) q[0];
rz(-0.61620617) q[1];
sx q[1];
rz(-2.1245978) q[1];
sx q[1];
rz(-2.9527051) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096029671) q[0];
sx q[0];
rz(-2.3052373) q[0];
sx q[0];
rz(1.0009074) q[0];
rz(-pi) q[1];
rz(-2.1470931) q[2];
sx q[2];
rz(-0.81953555) q[2];
sx q[2];
rz(2.799965) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9465883) q[1];
sx q[1];
rz(-2.431181) q[1];
sx q[1];
rz(-1.8844152) q[1];
rz(-pi) q[2];
rz(2.5152605) q[3];
sx q[3];
rz(-0.62343684) q[3];
sx q[3];
rz(-1.745519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6836493) q[2];
sx q[2];
rz(-1.8234437) q[2];
sx q[2];
rz(1.1434435) q[2];
rz(-2.5055366) q[3];
sx q[3];
rz(-3.0885185) q[3];
sx q[3];
rz(1.599357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77883333) q[0];
sx q[0];
rz(-2.9017359) q[0];
sx q[0];
rz(-1.973961) q[0];
rz(0.12655839) q[1];
sx q[1];
rz(-1.5153706) q[1];
sx q[1];
rz(-2.5680241) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1417057) q[0];
sx q[0];
rz(-0.83214271) q[0];
sx q[0];
rz(-3.0189199) q[0];
rz(2.5218256) q[2];
sx q[2];
rz(-2.7006442) q[2];
sx q[2];
rz(-0.82142219) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4497609) q[1];
sx q[1];
rz(-1.1402601) q[1];
sx q[1];
rz(3.0887763) q[1];
rz(-pi) q[2];
x q[2];
rz(0.6595936) q[3];
sx q[3];
rz(-1.9232188) q[3];
sx q[3];
rz(-0.11241684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6557189) q[2];
sx q[2];
rz(-1.8334917) q[2];
sx q[2];
rz(-1.6386848) q[2];
rz(1.8466628) q[3];
sx q[3];
rz(-0.27532268) q[3];
sx q[3];
rz(2.3143229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5824222) q[0];
sx q[0];
rz(-0.14635135) q[0];
sx q[0];
rz(2.225112) q[0];
rz(-2.9811663) q[1];
sx q[1];
rz(-1.8564936) q[1];
sx q[1];
rz(0.697335) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19577414) q[0];
sx q[0];
rz(-1.8911288) q[0];
sx q[0];
rz(2.8813502) q[0];
x q[1];
rz(-2.4795697) q[2];
sx q[2];
rz(-1.8141835) q[2];
sx q[2];
rz(3.0070729) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8749344) q[1];
sx q[1];
rz(-1.8968079) q[1];
sx q[1];
rz(-2.3951016) q[1];
rz(-pi) q[2];
rz(1.5847727) q[3];
sx q[3];
rz(-1.7009354) q[3];
sx q[3];
rz(-1.1443477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9852898) q[2];
sx q[2];
rz(-2.6142945) q[2];
sx q[2];
rz(-2.018833) q[2];
rz(-1.5491693) q[3];
sx q[3];
rz(-1.6808108) q[3];
sx q[3];
rz(0.39337015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8321946) q[0];
sx q[0];
rz(-0.10646146) q[0];
sx q[0];
rz(-1.1291809) q[0];
rz(0.87942266) q[1];
sx q[1];
rz(-1.8303454) q[1];
sx q[1];
rz(0.62854356) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.578772) q[0];
sx q[0];
rz(-1.6878753) q[0];
sx q[0];
rz(1.4463918) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4833335) q[2];
sx q[2];
rz(-1.1193917) q[2];
sx q[2];
rz(0.70027225) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.063911818) q[1];
sx q[1];
rz(-1.5305298) q[1];
sx q[1];
rz(-1.0185373) q[1];
rz(-pi) q[2];
rz(-0.83129779) q[3];
sx q[3];
rz(-1.6581088) q[3];
sx q[3];
rz(1.0333453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2200372) q[2];
sx q[2];
rz(-0.76801378) q[2];
sx q[2];
rz(-1.7680232) q[2];
rz(-2.8288362) q[3];
sx q[3];
rz(-1.0765358) q[3];
sx q[3];
rz(-0.086418644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3304928) q[0];
sx q[0];
rz(-2.436315) q[0];
sx q[0];
rz(0.16264859) q[0];
rz(1.9074408) q[1];
sx q[1];
rz(-1.1470497) q[1];
sx q[1];
rz(-2.2208234) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7013848) q[0];
sx q[0];
rz(-2.9281441) q[0];
sx q[0];
rz(2.0646854) q[0];
rz(-pi) q[1];
rz(0.55947907) q[2];
sx q[2];
rz(-2.1396416) q[2];
sx q[2];
rz(-0.43348962) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.28345767) q[1];
sx q[1];
rz(-2.3871536) q[1];
sx q[1];
rz(0.094875022) q[1];
x q[2];
rz(-0.1699888) q[3];
sx q[3];
rz(-2.3521479) q[3];
sx q[3];
rz(1.1657749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2938701) q[2];
sx q[2];
rz(-2.4203478) q[2];
sx q[2];
rz(-0.043962002) q[2];
rz(1.4219159) q[3];
sx q[3];
rz(-2.7094262) q[3];
sx q[3];
rz(-0.034984263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0111572) q[0];
sx q[0];
rz(-2.3373117) q[0];
sx q[0];
rz(-3.0479808) q[0];
rz(-1.0253819) q[1];
sx q[1];
rz(-1.5461812) q[1];
sx q[1];
rz(2.3536918) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073613361) q[0];
sx q[0];
rz(-0.96603051) q[0];
sx q[0];
rz(2.4725799) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56886832) q[2];
sx q[2];
rz(-1.7143357) q[2];
sx q[2];
rz(-1.6678866) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3026784) q[1];
sx q[1];
rz(-2.5414349) q[1];
sx q[1];
rz(0.072235302) q[1];
x q[2];
rz(1.4282088) q[3];
sx q[3];
rz(-1.6901921) q[3];
sx q[3];
rz(-2.6210331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0660144) q[2];
sx q[2];
rz(-2.5983073) q[2];
sx q[2];
rz(-0.96599609) q[2];
rz(-2.994359) q[3];
sx q[3];
rz(-1.4598264) q[3];
sx q[3];
rz(-2.537263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0886993) q[0];
sx q[0];
rz(-1.3857144) q[0];
sx q[0];
rz(1.2526441) q[0];
rz(-0.057627536) q[1];
sx q[1];
rz(-1.3725955) q[1];
sx q[1];
rz(2.7431814) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8114026) q[0];
sx q[0];
rz(-0.63985642) q[0];
sx q[0];
rz(-2.8485203) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6419027) q[2];
sx q[2];
rz(-1.4213287) q[2];
sx q[2];
rz(-0.88525822) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0115464) q[1];
sx q[1];
rz(-2.941644) q[1];
sx q[1];
rz(-1.5505095) q[1];
x q[2];
rz(1.4640973) q[3];
sx q[3];
rz(-1.0301642) q[3];
sx q[3];
rz(1.65834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5742699) q[2];
sx q[2];
rz(-1.8696573) q[2];
sx q[2];
rz(-2.2197913) q[2];
rz(0.70720339) q[3];
sx q[3];
rz(-1.0450398) q[3];
sx q[3];
rz(-0.31064492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9486174) q[0];
sx q[0];
rz(-3.0204168) q[0];
sx q[0];
rz(-0.90712547) q[0];
rz(-0.47997296) q[1];
sx q[1];
rz(-2.0860806) q[1];
sx q[1];
rz(2.5352246) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84068882) q[0];
sx q[0];
rz(-0.32561526) q[0];
sx q[0];
rz(-0.31890829) q[0];
rz(2.3263518) q[2];
sx q[2];
rz(-2.4324907) q[2];
sx q[2];
rz(-0.39295024) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.94928259) q[1];
sx q[1];
rz(-2.3282166) q[1];
sx q[1];
rz(1.1622914) q[1];
rz(-pi) q[2];
rz(0.65878792) q[3];
sx q[3];
rz(-1.4843936) q[3];
sx q[3];
rz(0.98443809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.85764641) q[2];
sx q[2];
rz(-1.5479156) q[2];
sx q[2];
rz(-0.66961163) q[2];
rz(-2.5745463) q[3];
sx q[3];
rz(-1.0457057) q[3];
sx q[3];
rz(-1.4001575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46212101) q[0];
sx q[0];
rz(-1.7993878) q[0];
sx q[0];
rz(-0.8959499) q[0];
rz(3.1013007) q[1];
sx q[1];
rz(-2.1998684) q[1];
sx q[1];
rz(-1.5641854) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58916726) q[0];
sx q[0];
rz(-1.5755113) q[0];
sx q[0];
rz(1.5472359) q[0];
x q[1];
rz(-2.0032004) q[2];
sx q[2];
rz(-0.93195019) q[2];
sx q[2];
rz(2.2496719) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2103232) q[1];
sx q[1];
rz(-2.0788553) q[1];
sx q[1];
rz(2.6711051) q[1];
rz(-pi) q[2];
rz(0.96503105) q[3];
sx q[3];
rz(-1.8739227) q[3];
sx q[3];
rz(-2.3245963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3080052) q[2];
sx q[2];
rz(-2.8218125) q[2];
sx q[2];
rz(-1.3099111) q[2];
rz(-1.9561907) q[3];
sx q[3];
rz(-2.2535321) q[3];
sx q[3];
rz(1.5230702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1401405) q[0];
sx q[0];
rz(-1.3012713) q[0];
sx q[0];
rz(-0.95300994) q[0];
rz(3.0630655) q[1];
sx q[1];
rz(-2.0546866) q[1];
sx q[1];
rz(2.3436117) q[1];
rz(-0.8962411) q[2];
sx q[2];
rz(-1.3660539) q[2];
sx q[2];
rz(2.9851253) q[2];
rz(-2.9361712) q[3];
sx q[3];
rz(-1.724549) q[3];
sx q[3];
rz(-2.1461501) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
