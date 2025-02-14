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
rz(-1.4650605) q[0];
sx q[0];
rz(-0.2413916) q[0];
sx q[0];
rz(-3.0093332) q[0];
rz(-0.013068696) q[1];
sx q[1];
rz(-0.71813923) q[1];
sx q[1];
rz(-0.018996039) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26532224) q[0];
sx q[0];
rz(-0.92060584) q[0];
sx q[0];
rz(0.27133915) q[0];
rz(-pi) q[1];
rz(-1.1055345) q[2];
sx q[2];
rz(-0.20814116) q[2];
sx q[2];
rz(0.15424745) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2680696) q[1];
sx q[1];
rz(-1.8933081) q[1];
sx q[1];
rz(-1.3871865) q[1];
x q[2];
rz(-2.1601342) q[3];
sx q[3];
rz(-2.8638726) q[3];
sx q[3];
rz(-0.026175682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.12048177) q[2];
sx q[2];
rz(-2.3349473) q[2];
sx q[2];
rz(-0.79180229) q[2];
rz(1.0857438) q[3];
sx q[3];
rz(-1.9224242) q[3];
sx q[3];
rz(-2.048548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.87078142) q[0];
sx q[0];
rz(-2.9830611) q[0];
sx q[0];
rz(-2.812401) q[0];
rz(1.6393433) q[1];
sx q[1];
rz(-0.90198016) q[1];
sx q[1];
rz(0.42218581) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6471408) q[0];
sx q[0];
rz(-1.1472771) q[0];
sx q[0];
rz(-0.5943055) q[0];
rz(-pi) q[1];
rz(0.54791003) q[2];
sx q[2];
rz(-1.413125) q[2];
sx q[2];
rz(1.6239177) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.364155) q[1];
sx q[1];
rz(-1.5244487) q[1];
sx q[1];
rz(-0.32407659) q[1];
x q[2];
rz(-3.0358697) q[3];
sx q[3];
rz(-0.7751677) q[3];
sx q[3];
rz(-0.4251484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.72055703) q[2];
sx q[2];
rz(-1.9562419) q[2];
sx q[2];
rz(-1.1400918) q[2];
rz(-3.1371878) q[3];
sx q[3];
rz(-1.5811698) q[3];
sx q[3];
rz(0.13979039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7797101) q[0];
sx q[0];
rz(-3.0411868) q[0];
sx q[0];
rz(0.34565872) q[0];
rz(2.0480305) q[1];
sx q[1];
rz(-0.82229096) q[1];
sx q[1];
rz(-2.0872769) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38906583) q[0];
sx q[0];
rz(-1.1841762) q[0];
sx q[0];
rz(1.9324612) q[0];
x q[1];
rz(0.81297154) q[2];
sx q[2];
rz(-1.6457498) q[2];
sx q[2];
rz(3.093442) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.24040996) q[1];
sx q[1];
rz(-0.4571656) q[1];
sx q[1];
rz(0.14480503) q[1];
rz(-pi) q[2];
rz(3.1183447) q[3];
sx q[3];
rz(-2.23684) q[3];
sx q[3];
rz(2.9262528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.67636079) q[2];
sx q[2];
rz(-2.5939442) q[2];
sx q[2];
rz(-2.5993627) q[2];
rz(-0.17255653) q[3];
sx q[3];
rz(-1.4901284) q[3];
sx q[3];
rz(2.0358613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3365823) q[0];
sx q[0];
rz(-1.35291) q[0];
sx q[0];
rz(-1.8804869) q[0];
rz(1.1359967) q[1];
sx q[1];
rz(-2.0295862) q[1];
sx q[1];
rz(0.13066185) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23993429) q[0];
sx q[0];
rz(-1.8404418) q[0];
sx q[0];
rz(1.1872227) q[0];
rz(2.7515) q[2];
sx q[2];
rz(-0.54020262) q[2];
sx q[2];
rz(-2.982123) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9798292) q[1];
sx q[1];
rz(-1.5391304) q[1];
sx q[1];
rz(2.0791441) q[1];
rz(-0.66417905) q[3];
sx q[3];
rz(-1.8669493) q[3];
sx q[3];
rz(-0.68063762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7589492) q[2];
sx q[2];
rz(-0.073315695) q[2];
sx q[2];
rz(-2.9200413) q[2];
rz(-2.4394636) q[3];
sx q[3];
rz(-1.5626855) q[3];
sx q[3];
rz(2.4041972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.57109433) q[0];
sx q[0];
rz(-1.8966738) q[0];
sx q[0];
rz(3.0076497) q[0];
rz(-0.36930034) q[1];
sx q[1];
rz(-1.1993273) q[1];
sx q[1];
rz(1.7514924) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54467359) q[0];
sx q[0];
rz(-1.8159165) q[0];
sx q[0];
rz(-1.0463723) q[0];
rz(-pi) q[1];
rz(0.26525396) q[2];
sx q[2];
rz(-1.2758288) q[2];
sx q[2];
rz(-0.8910999) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.578957) q[1];
sx q[1];
rz(-1.0308415) q[1];
sx q[1];
rz(-1.8946429) q[1];
rz(-0.36365328) q[3];
sx q[3];
rz(-1.8580274) q[3];
sx q[3];
rz(-1.2032697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.2964581) q[2];
sx q[2];
rz(-0.90201169) q[2];
sx q[2];
rz(1.585539) q[2];
rz(0.29160094) q[3];
sx q[3];
rz(-2.6652938) q[3];
sx q[3];
rz(0.27879032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6911102) q[0];
sx q[0];
rz(-2.1463558) q[0];
sx q[0];
rz(2.6158748) q[0];
rz(-2.763343) q[1];
sx q[1];
rz(-1.6098166) q[1];
sx q[1];
rz(2.8674616) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0167086) q[0];
sx q[0];
rz(-1.5450461) q[0];
sx q[0];
rz(0.59016778) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1172764) q[2];
sx q[2];
rz(-1.9651061) q[2];
sx q[2];
rz(2.0298983) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.25933274) q[1];
sx q[1];
rz(-1.3162398) q[1];
sx q[1];
rz(-1.1544636) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.317362) q[3];
sx q[3];
rz(-1.715664) q[3];
sx q[3];
rz(-1.0235909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5963001) q[2];
sx q[2];
rz(-1.9175074) q[2];
sx q[2];
rz(2.8833585) q[2];
rz(-1.5457414) q[3];
sx q[3];
rz(-1.7786547) q[3];
sx q[3];
rz(-0.405092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79893583) q[0];
sx q[0];
rz(-2.6234143) q[0];
sx q[0];
rz(-0.10547353) q[0];
rz(-2.8666829) q[1];
sx q[1];
rz(-1.4242438) q[1];
sx q[1];
rz(-2.8555433) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8409922) q[0];
sx q[0];
rz(-1.6874041) q[0];
sx q[0];
rz(-2.9564315) q[0];
rz(-pi) q[1];
rz(1.8931383) q[2];
sx q[2];
rz(-2.3791056) q[2];
sx q[2];
rz(0.74021268) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6374928) q[1];
sx q[1];
rz(-1.6352202) q[1];
sx q[1];
rz(1.9570051) q[1];
rz(2.5887542) q[3];
sx q[3];
rz(-1.1854608) q[3];
sx q[3];
rz(-0.60147731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1470571) q[2];
sx q[2];
rz(-2.4592082) q[2];
sx q[2];
rz(0.4500173) q[2];
rz(-0.58722812) q[3];
sx q[3];
rz(-1.3415895) q[3];
sx q[3];
rz(-1.8915061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8477972) q[0];
sx q[0];
rz(-1.6814517) q[0];
sx q[0];
rz(-1.1867123) q[0];
rz(2.8222491) q[1];
sx q[1];
rz(-2.1198544) q[1];
sx q[1];
rz(0.26225463) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38459331) q[0];
sx q[0];
rz(-0.5890863) q[0];
sx q[0];
rz(0.37269816) q[0];
rz(0.13599254) q[2];
sx q[2];
rz(-1.9806702) q[2];
sx q[2];
rz(2.6213344) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9567392) q[1];
sx q[1];
rz(-1.0295086) q[1];
sx q[1];
rz(-0.095603099) q[1];
rz(0.45164177) q[3];
sx q[3];
rz(-1.9596976) q[3];
sx q[3];
rz(1.1448142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.35683262) q[2];
sx q[2];
rz(-0.46669745) q[2];
sx q[2];
rz(-2.812815) q[2];
rz(1.6262936) q[3];
sx q[3];
rz(-1.7833775) q[3];
sx q[3];
rz(2.7330107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3560155) q[0];
sx q[0];
rz(-2.6308036) q[0];
sx q[0];
rz(0.27895862) q[0];
rz(0.51697671) q[1];
sx q[1];
rz(-0.37716436) q[1];
sx q[1];
rz(1.4252211) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12869975) q[0];
sx q[0];
rz(-2.1919247) q[0];
sx q[0];
rz(-1.6245882) q[0];
rz(-2.4819229) q[2];
sx q[2];
rz(-2.2222509) q[2];
sx q[2];
rz(1.720681) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9659075) q[1];
sx q[1];
rz(-2.0669427) q[1];
sx q[1];
rz(1.7860852) q[1];
rz(-pi) q[2];
rz(-2.2278085) q[3];
sx q[3];
rz(-1.3385286) q[3];
sx q[3];
rz(-2.6871347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4031389) q[2];
sx q[2];
rz(-2.8454915) q[2];
sx q[2];
rz(2.9816755) q[2];
rz(-0.81965172) q[3];
sx q[3];
rz(-1.22217) q[3];
sx q[3];
rz(2.405449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4495471) q[0];
sx q[0];
rz(-1.8501546) q[0];
sx q[0];
rz(-2.8458169) q[0];
rz(-2.6462818) q[1];
sx q[1];
rz(-0.43423978) q[1];
sx q[1];
rz(1.6326509) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42351022) q[0];
sx q[0];
rz(-2.5959466) q[0];
sx q[0];
rz(-0.83828022) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7118218) q[2];
sx q[2];
rz(-2.786028) q[2];
sx q[2];
rz(2.9586458) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.26824046) q[1];
sx q[1];
rz(-0.11246364) q[1];
sx q[1];
rz(-2.4960773) q[1];
rz(-pi) q[2];
rz(-1.1811352) q[3];
sx q[3];
rz(-1.1077322) q[3];
sx q[3];
rz(0.72805007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.32790023) q[2];
sx q[2];
rz(-1.6715965) q[2];
sx q[2];
rz(0.63391614) q[2];
rz(3.0564803) q[3];
sx q[3];
rz(-1.8952993) q[3];
sx q[3];
rz(2.3688721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2278628) q[0];
sx q[0];
rz(-1.5416523) q[0];
sx q[0];
rz(-0.1097485) q[0];
rz(1.2211424) q[1];
sx q[1];
rz(-2.2962062) q[1];
sx q[1];
rz(-2.5795945) q[1];
rz(1.1176422) q[2];
sx q[2];
rz(-0.48095555) q[2];
sx q[2];
rz(-1.7775735) q[2];
rz(3.1172995) q[3];
sx q[3];
rz(-1.3391936) q[3];
sx q[3];
rz(-1.5486123) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
