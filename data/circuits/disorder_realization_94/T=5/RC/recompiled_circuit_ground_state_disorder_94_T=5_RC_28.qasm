OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2151467) q[0];
sx q[0];
rz(-2.4617221) q[0];
sx q[0];
rz(-3.0409066) q[0];
rz(0.45473948) q[1];
sx q[1];
rz(-0.68128959) q[1];
sx q[1];
rz(-1.0256306) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89698638) q[0];
sx q[0];
rz(-1.7660487) q[0];
sx q[0];
rz(0.050873916) q[0];
x q[1];
rz(2.8335613) q[2];
sx q[2];
rz(-1.6197761) q[2];
sx q[2];
rz(2.2394951) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.611022) q[1];
sx q[1];
rz(-0.82029063) q[1];
sx q[1];
rz(-1.6802701) q[1];
rz(-pi) q[2];
rz(2.1061534) q[3];
sx q[3];
rz(-0.97781721) q[3];
sx q[3];
rz(-1.6361332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2166298) q[2];
sx q[2];
rz(-1.7731885) q[2];
sx q[2];
rz(-1.487757) q[2];
rz(1.1612859) q[3];
sx q[3];
rz(-2.1372644) q[3];
sx q[3];
rz(1.0950834) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0733114) q[0];
sx q[0];
rz(-0.40638766) q[0];
sx q[0];
rz(2.708129) q[0];
rz(-2.0794226) q[1];
sx q[1];
rz(-1.559161) q[1];
sx q[1];
rz(-0.72057048) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67340785) q[0];
sx q[0];
rz(-1.5338729) q[0];
sx q[0];
rz(1.6975178) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8938204) q[2];
sx q[2];
rz(-1.0508014) q[2];
sx q[2];
rz(-1.8682075) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7097672) q[1];
sx q[1];
rz(-2.7569237) q[1];
sx q[1];
rz(2.7363214) q[1];
rz(1.0508762) q[3];
sx q[3];
rz(-2.5934016) q[3];
sx q[3];
rz(0.91778008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7945107) q[2];
sx q[2];
rz(-1.2257267) q[2];
sx q[2];
rz(2.1573055) q[2];
rz(0.85683626) q[3];
sx q[3];
rz(-0.58647668) q[3];
sx q[3];
rz(-1.8872567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1089351) q[0];
sx q[0];
rz(-0.28626838) q[0];
sx q[0];
rz(-0.74288595) q[0];
rz(-3.1318829) q[1];
sx q[1];
rz(-1.5747036) q[1];
sx q[1];
rz(2.8579393) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46732908) q[0];
sx q[0];
rz(-2.0785851) q[0];
sx q[0];
rz(-1.171815) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.1998603) q[2];
sx q[2];
rz(-1.4132573) q[2];
sx q[2];
rz(0.021856088) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.7934809) q[1];
sx q[1];
rz(-1.8758131) q[1];
sx q[1];
rz(-2.9654337) q[1];
rz(-pi) q[2];
rz(-3.1183692) q[3];
sx q[3];
rz(-1.4693714) q[3];
sx q[3];
rz(-1.2161127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.11188406) q[2];
sx q[2];
rz(-1.6969029) q[2];
sx q[2];
rz(0.65579826) q[2];
rz(-0.2119952) q[3];
sx q[3];
rz(-2.5097804) q[3];
sx q[3];
rz(-1.9688212) q[3];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028932171) q[0];
sx q[0];
rz(-2.6472968) q[0];
sx q[0];
rz(-1.5546881) q[0];
rz(-0.2298062) q[1];
sx q[1];
rz(-2.4202012) q[1];
sx q[1];
rz(-1.8386286) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7763846) q[0];
sx q[0];
rz(-1.0321277) q[0];
sx q[0];
rz(-2.0950277) q[0];
rz(0.26891687) q[2];
sx q[2];
rz(-1.8559716) q[2];
sx q[2];
rz(-2.0320924) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.069217056) q[1];
sx q[1];
rz(-2.0942959) q[1];
sx q[1];
rz(-1.545946) q[1];
rz(-1.0082208) q[3];
sx q[3];
rz(-1.8726714) q[3];
sx q[3];
rz(-2.9396597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.135123) q[2];
sx q[2];
rz(-2.3389356) q[2];
sx q[2];
rz(2.6378677) q[2];
rz(-1.6040365) q[3];
sx q[3];
rz(-1.9683014) q[3];
sx q[3];
rz(1.5639719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2876494) q[0];
sx q[0];
rz(-2.8048153) q[0];
sx q[0];
rz(2.7119998) q[0];
rz(-2.8751539) q[1];
sx q[1];
rz(-0.8100422) q[1];
sx q[1];
rz(2.0727167) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3392068) q[0];
sx q[0];
rz(-0.9046208) q[0];
sx q[0];
rz(1.0399014) q[0];
rz(-pi) q[1];
rz(-1.2524302) q[2];
sx q[2];
rz(-2.3471544) q[2];
sx q[2];
rz(2.3311262) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7366745) q[1];
sx q[1];
rz(-2.6231986) q[1];
sx q[1];
rz(0.53804086) q[1];
x q[2];
rz(-0.74917082) q[3];
sx q[3];
rz(-0.35532829) q[3];
sx q[3];
rz(-0.65566777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2499007) q[2];
sx q[2];
rz(-0.72600681) q[2];
sx q[2];
rz(0.1740087) q[2];
rz(-1.2627259) q[3];
sx q[3];
rz(-1.5401253) q[3];
sx q[3];
rz(-1.5549829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0083017666) q[0];
sx q[0];
rz(-1.0187848) q[0];
sx q[0];
rz(0.89023501) q[0];
rz(1.4324191) q[1];
sx q[1];
rz(-1.018367) q[1];
sx q[1];
rz(0.11996809) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.650506) q[0];
sx q[0];
rz(-1.0179369) q[0];
sx q[0];
rz(-1.4006459) q[0];
rz(-2.2416297) q[2];
sx q[2];
rz(-2.6113178) q[2];
sx q[2];
rz(1.8535623) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.97291291) q[1];
sx q[1];
rz(-0.69818234) q[1];
sx q[1];
rz(0.38777988) q[1];
rz(-2.6011068) q[3];
sx q[3];
rz(-1.4048164) q[3];
sx q[3];
rz(-2.9946208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1098108) q[2];
sx q[2];
rz(-1.6620771) q[2];
sx q[2];
rz(0.39153448) q[2];
rz(-1.3854965) q[3];
sx q[3];
rz(-2.9031495) q[3];
sx q[3];
rz(2.5029206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2077797) q[0];
sx q[0];
rz(-0.55140984) q[0];
sx q[0];
rz(-1.3054003) q[0];
rz(0.50453672) q[1];
sx q[1];
rz(-1.7326109) q[1];
sx q[1];
rz(-0.22444589) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69687122) q[0];
sx q[0];
rz(-1.5550647) q[0];
sx q[0];
rz(-0.016642668) q[0];
rz(-pi) q[1];
rz(1.7845909) q[2];
sx q[2];
rz(-2.0540621) q[2];
sx q[2];
rz(2.4501806) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.689381) q[1];
sx q[1];
rz(-1.2465917) q[1];
sx q[1];
rz(2.8597742) q[1];
rz(-pi) q[2];
rz(-1.2013632) q[3];
sx q[3];
rz(-1.2681586) q[3];
sx q[3];
rz(1.7517881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9597783) q[2];
sx q[2];
rz(-2.0373127) q[2];
sx q[2];
rz(-0.9474729) q[2];
rz(3.0090561) q[3];
sx q[3];
rz(-0.65968502) q[3];
sx q[3];
rz(0.98602492) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2333616) q[0];
sx q[0];
rz(-2.7955604) q[0];
sx q[0];
rz(-2.76873) q[0];
rz(2.8064959) q[1];
sx q[1];
rz(-1.9357888) q[1];
sx q[1];
rz(-2.0705409) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11075739) q[0];
sx q[0];
rz(-1.6399151) q[0];
sx q[0];
rz(-1.5821619) q[0];
x q[1];
rz(-0.55166371) q[2];
sx q[2];
rz(-1.3886189) q[2];
sx q[2];
rz(-2.7523486) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.95360699) q[1];
sx q[1];
rz(-1.803557) q[1];
sx q[1];
rz(-2.8899419) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2921807) q[3];
sx q[3];
rz(-2.8077586) q[3];
sx q[3];
rz(-2.1899471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4272473) q[2];
sx q[2];
rz(-1.599396) q[2];
sx q[2];
rz(2.074312) q[2];
rz(0.37211564) q[3];
sx q[3];
rz(-0.99388638) q[3];
sx q[3];
rz(-2.9414162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9099971) q[0];
sx q[0];
rz(-0.59760004) q[0];
sx q[0];
rz(-1.4071314) q[0];
rz(2.097997) q[1];
sx q[1];
rz(-0.93799543) q[1];
sx q[1];
rz(-0.49120206) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.585116) q[0];
sx q[0];
rz(-2.399754) q[0];
sx q[0];
rz(2.0061352) q[0];
rz(0.78989002) q[2];
sx q[2];
rz(-0.23755506) q[2];
sx q[2];
rz(1.6611163) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.421153) q[1];
sx q[1];
rz(-0.3657839) q[1];
sx q[1];
rz(0.37558742) q[1];
x q[2];
rz(1.7513428) q[3];
sx q[3];
rz(-1.9420764) q[3];
sx q[3];
rz(-0.89902395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6970984) q[2];
sx q[2];
rz(-1.5105057) q[2];
sx q[2];
rz(2.0972882) q[2];
rz(-0.77110428) q[3];
sx q[3];
rz(-1.5325357) q[3];
sx q[3];
rz(0.3573629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
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
rz(-0.20030178) q[0];
sx q[0];
rz(-1.8481978) q[0];
sx q[0];
rz(-1.6218761) q[0];
rz(-1.8408076) q[1];
sx q[1];
rz(-1.3011353) q[1];
sx q[1];
rz(0.69449743) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3169169) q[0];
sx q[0];
rz(-1.7986713) q[0];
sx q[0];
rz(0.014046305) q[0];
rz(-pi) q[1];
rz(1.7037665) q[2];
sx q[2];
rz(-2.0914234) q[2];
sx q[2];
rz(-2.7085053) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4150691) q[1];
sx q[1];
rz(-1.5818672) q[1];
sx q[1];
rz(1.2367718) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1841838) q[3];
sx q[3];
rz(-2.2848633) q[3];
sx q[3];
rz(2.3165645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0684315) q[2];
sx q[2];
rz(-1.3621) q[2];
sx q[2];
rz(2.6279602) q[2];
rz(2.8237776) q[3];
sx q[3];
rz(-1.5376667) q[3];
sx q[3];
rz(-1.0191127) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1401055) q[0];
sx q[0];
rz(-1.5167863) q[0];
sx q[0];
rz(2.232502) q[0];
rz(-1.2263251) q[1];
sx q[1];
rz(-1.6188123) q[1];
sx q[1];
rz(0.34368044) q[1];
rz(2.0991391) q[2];
sx q[2];
rz(-2.1975542) q[2];
sx q[2];
rz(2.2007813) q[2];
rz(1.574434) q[3];
sx q[3];
rz(-0.52593492) q[3];
sx q[3];
rz(1.8686663) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
