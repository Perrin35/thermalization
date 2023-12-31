OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2919579) q[0];
sx q[0];
rz(-2.7014974) q[0];
sx q[0];
rz(-0.13719288) q[0];
rz(-1.7358915) q[1];
sx q[1];
rz(-1.403221) q[1];
sx q[1];
rz(-0.52991968) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7512902) q[0];
sx q[0];
rz(-1.6781224) q[0];
sx q[0];
rz(-1.1885378) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9723188) q[2];
sx q[2];
rz(-1.8895154) q[2];
sx q[2];
rz(-2.4674494) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4094028) q[1];
sx q[1];
rz(-2.4362872) q[1];
sx q[1];
rz(0.7440872) q[1];
rz(-pi) q[2];
rz(2.3832541) q[3];
sx q[3];
rz(-1.7539277) q[3];
sx q[3];
rz(-1.5116215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4522176) q[2];
sx q[2];
rz(-1.3000501) q[2];
sx q[2];
rz(-2.8049862) q[2];
rz(-1.5161139) q[3];
sx q[3];
rz(-0.55364004) q[3];
sx q[3];
rz(-1.6158993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34823725) q[0];
sx q[0];
rz(-1.1084778) q[0];
sx q[0];
rz(-3.120378) q[0];
rz(-1.1938098) q[1];
sx q[1];
rz(-2.1021011) q[1];
sx q[1];
rz(-2.3056727) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0854557) q[0];
sx q[0];
rz(-1.5888927) q[0];
sx q[0];
rz(1.7096814) q[0];
rz(2.1199273) q[2];
sx q[2];
rz(-1.9252535) q[2];
sx q[2];
rz(-2.3159388) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6902496) q[1];
sx q[1];
rz(-0.93938821) q[1];
sx q[1];
rz(2.8020225) q[1];
rz(0.73987506) q[3];
sx q[3];
rz(-1.8679108) q[3];
sx q[3];
rz(-1.5418996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2772284) q[2];
sx q[2];
rz(-1.9936864) q[2];
sx q[2];
rz(-1.345984) q[2];
rz(-0.35955444) q[3];
sx q[3];
rz(-0.94272009) q[3];
sx q[3];
rz(-0.49697044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7132752) q[0];
sx q[0];
rz(-2.064216) q[0];
sx q[0];
rz(1.0536449) q[0];
rz(-1.2288278) q[1];
sx q[1];
rz(-1.5412953) q[1];
sx q[1];
rz(-2.704481) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8489704) q[0];
sx q[0];
rz(-1.4842352) q[0];
sx q[0];
rz(1.2785046) q[0];
x q[1];
rz(-2.897981) q[2];
sx q[2];
rz(-1.9677475) q[2];
sx q[2];
rz(-1.5543907) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6091842) q[1];
sx q[1];
rz(-0.44645616) q[1];
sx q[1];
rz(-0.32534728) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0655572) q[3];
sx q[3];
rz(-1.6238188) q[3];
sx q[3];
rz(-0.094129063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.019471021) q[2];
sx q[2];
rz(-0.78142587) q[2];
sx q[2];
rz(-2.1195228) q[2];
rz(-1.9034889) q[3];
sx q[3];
rz(-0.3823897) q[3];
sx q[3];
rz(-0.4195655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5220752) q[0];
sx q[0];
rz(-1.2457122) q[0];
sx q[0];
rz(-0.98130256) q[0];
rz(3.006382) q[1];
sx q[1];
rz(-1.0842666) q[1];
sx q[1];
rz(-2.9503126) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0606196) q[0];
sx q[0];
rz(-1.4220211) q[0];
sx q[0];
rz(0.087555126) q[0];
rz(0.27387597) q[2];
sx q[2];
rz(-2.2857776) q[2];
sx q[2];
rz(2.3194734) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3595708) q[1];
sx q[1];
rz(-1.6732209) q[1];
sx q[1];
rz(-0.77918474) q[1];
rz(-pi) q[2];
rz(2.5883834) q[3];
sx q[3];
rz(-0.60713967) q[3];
sx q[3];
rz(0.34724423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4613351) q[2];
sx q[2];
rz(-0.98538435) q[2];
sx q[2];
rz(-2.130924) q[2];
rz(-2.3800395) q[3];
sx q[3];
rz(-1.1798309) q[3];
sx q[3];
rz(0.23553577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33870944) q[0];
sx q[0];
rz(-0.25512472) q[0];
sx q[0];
rz(-2.5849735) q[0];
rz(3.026475) q[1];
sx q[1];
rz(-1.8042253) q[1];
sx q[1];
rz(-2.1690878) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92161979) q[0];
sx q[0];
rz(-0.12244206) q[0];
sx q[0];
rz(0.58453154) q[0];
rz(-2.3210578) q[2];
sx q[2];
rz(-1.3247326) q[2];
sx q[2];
rz(2.8503502) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.79341187) q[1];
sx q[1];
rz(-0.23985292) q[1];
sx q[1];
rz(-2.1972149) q[1];
rz(3.1049018) q[3];
sx q[3];
rz(-2.1864236) q[3];
sx q[3];
rz(2.7021367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.82289034) q[2];
sx q[2];
rz(-1.0914785) q[2];
sx q[2];
rz(-1.548432) q[2];
rz(1.7758153) q[3];
sx q[3];
rz(-0.32309353) q[3];
sx q[3];
rz(2.1877066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4218629) q[0];
sx q[0];
rz(-1.8780163) q[0];
sx q[0];
rz(1.4259889) q[0];
rz(-1.0643719) q[1];
sx q[1];
rz(-2.1247037) q[1];
sx q[1];
rz(2.7672966) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4842589) q[0];
sx q[0];
rz(-1.1614292) q[0];
sx q[0];
rz(-1.9166458) q[0];
rz(-pi) q[1];
rz(-2.6815368) q[2];
sx q[2];
rz(-0.54121491) q[2];
sx q[2];
rz(0.40749007) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6090138) q[1];
sx q[1];
rz(-2.5206869) q[1];
sx q[1];
rz(0.46355526) q[1];
x q[2];
rz(-0.66773325) q[3];
sx q[3];
rz(-2.1453834) q[3];
sx q[3];
rz(-2.7001065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8950243) q[2];
sx q[2];
rz(-0.66528577) q[2];
sx q[2];
rz(-0.95823112) q[2];
rz(0.22917497) q[3];
sx q[3];
rz(-1.6848247) q[3];
sx q[3];
rz(-0.62098256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(0.36528698) q[0];
sx q[0];
rz(-1.1927274) q[0];
sx q[0];
rz(2.2348485) q[0];
rz(2.0523741) q[1];
sx q[1];
rz(-1.6420495) q[1];
sx q[1];
rz(1.8315171) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.487264) q[0];
sx q[0];
rz(-1.9582821) q[0];
sx q[0];
rz(1.5154454) q[0];
rz(-1.2398948) q[2];
sx q[2];
rz(-2.2563997) q[2];
sx q[2];
rz(2.9943525) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8430427) q[1];
sx q[1];
rz(-1.7610465) q[1];
sx q[1];
rz(-2.802554) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0974) q[3];
sx q[3];
rz(-0.5558388) q[3];
sx q[3];
rz(-1.0122055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.92581302) q[2];
sx q[2];
rz(-0.44162193) q[2];
sx q[2];
rz(-1.6581992) q[2];
rz(0.27967927) q[3];
sx q[3];
rz(-2.1488583) q[3];
sx q[3];
rz(-0.057597615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16335547) q[0];
sx q[0];
rz(-1.6219448) q[0];
sx q[0];
rz(-2.9220007) q[0];
rz(-2.638468) q[1];
sx q[1];
rz(-0.88880912) q[1];
sx q[1];
rz(-0.84987744) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6760611) q[0];
sx q[0];
rz(-1.7277576) q[0];
sx q[0];
rz(1.7992875) q[0];
rz(-pi) q[1];
rz(-1.6860028) q[2];
sx q[2];
rz(-0.68636471) q[2];
sx q[2];
rz(-0.79730588) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1666607) q[1];
sx q[1];
rz(-2.3798675) q[1];
sx q[1];
rz(1.6208266) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0073118) q[3];
sx q[3];
rz(-2.8562162) q[3];
sx q[3];
rz(-0.34261045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5510817) q[2];
sx q[2];
rz(-1.8293646) q[2];
sx q[2];
rz(1.760651) q[2];
rz(-2.3855709) q[3];
sx q[3];
rz(-0.20320007) q[3];
sx q[3];
rz(0.35593629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6324156) q[0];
sx q[0];
rz(-2.248705) q[0];
sx q[0];
rz(0.40503043) q[0];
rz(-0.45267725) q[1];
sx q[1];
rz(-0.98222268) q[1];
sx q[1];
rz(1.2776432) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1179827) q[0];
sx q[0];
rz(-1.3408459) q[0];
sx q[0];
rz(-2.635637) q[0];
rz(2.5759376) q[2];
sx q[2];
rz(-0.7753765) q[2];
sx q[2];
rz(-2.2765809) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8511508) q[1];
sx q[1];
rz(-2.6944707) q[1];
sx q[1];
rz(1.4852344) q[1];
x q[2];
rz(2.7550335) q[3];
sx q[3];
rz(-0.87214008) q[3];
sx q[3];
rz(-0.78702918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3383125) q[2];
sx q[2];
rz(-0.40643224) q[2];
sx q[2];
rz(2.7837616) q[2];
rz(1.4194277) q[3];
sx q[3];
rz(-1.8678886) q[3];
sx q[3];
rz(2.0675802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(-2.2232067) q[0];
sx q[0];
rz(-3.0637488) q[0];
sx q[0];
rz(0.11225587) q[0];
rz(-2.2414801) q[1];
sx q[1];
rz(-1.0670412) q[1];
sx q[1];
rz(2.9311438) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1505517) q[0];
sx q[0];
rz(-1.7518839) q[0];
sx q[0];
rz(-2.6291356) q[0];
rz(0.16410447) q[2];
sx q[2];
rz(-1.9774984) q[2];
sx q[2];
rz(-1.1526398) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.12293832) q[1];
sx q[1];
rz(-1.2631053) q[1];
sx q[1];
rz(-1.3535181) q[1];
rz(-pi) q[2];
rz(2.3841303) q[3];
sx q[3];
rz(-2.8031581) q[3];
sx q[3];
rz(0.38871845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9615053) q[2];
sx q[2];
rz(-2.4752361) q[2];
sx q[2];
rz(-1.5562742) q[2];
rz(1.8680343) q[3];
sx q[3];
rz(-0.62265101) q[3];
sx q[3];
rz(-0.20726985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56959854) q[0];
sx q[0];
rz(-0.8710237) q[0];
sx q[0];
rz(-1.3652753) q[0];
rz(-2.3251484) q[1];
sx q[1];
rz(-1.888231) q[1];
sx q[1];
rz(2.9838557) q[1];
rz(1.1007166) q[2];
sx q[2];
rz(-1.7190949) q[2];
sx q[2];
rz(0.20863056) q[2];
rz(-2.7995085) q[3];
sx q[3];
rz(-0.87089201) q[3];
sx q[3];
rz(1.0637829) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
