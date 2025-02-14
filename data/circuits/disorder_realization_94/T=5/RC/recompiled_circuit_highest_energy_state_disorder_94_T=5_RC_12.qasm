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
rz(0.87297451) q[0];
sx q[0];
rz(-0.31384808) q[0];
sx q[0];
rz(2.0883972) q[0];
rz(-1.516951) q[1];
sx q[1];
rz(-2.8487974) q[1];
sx q[1];
rz(2.204978) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7335986) q[0];
sx q[0];
rz(-0.50845611) q[0];
sx q[0];
rz(-1.5336799) q[0];
rz(-pi) q[1];
rz(0.99586113) q[2];
sx q[2];
rz(-2.323192) q[2];
sx q[2];
rz(-0.40413293) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.577212) q[1];
sx q[1];
rz(-2.138277) q[1];
sx q[1];
rz(1.903731) q[1];
rz(-pi) q[2];
rz(-2.7500912) q[3];
sx q[3];
rz(-1.748852) q[3];
sx q[3];
rz(-0.72339971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8006353) q[2];
sx q[2];
rz(-2.1850047) q[2];
sx q[2];
rz(-2.0802278) q[2];
rz(0.64570767) q[3];
sx q[3];
rz(-2.784745) q[3];
sx q[3];
rz(2.1319353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89479947) q[0];
sx q[0];
rz(-0.49870393) q[0];
sx q[0];
rz(0.45951581) q[0];
rz(-1.3872321) q[1];
sx q[1];
rz(-0.51001716) q[1];
sx q[1];
rz(2.0769108) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6094006) q[0];
sx q[0];
rz(-1.4123035) q[0];
sx q[0];
rz(-2.81937) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26724813) q[2];
sx q[2];
rz(-1.5540872) q[2];
sx q[2];
rz(0.087053336) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3749417) q[1];
sx q[1];
rz(-1.8901575) q[1];
sx q[1];
rz(1.6691895) q[1];
x q[2];
rz(1.2336938) q[3];
sx q[3];
rz(-1.8300984) q[3];
sx q[3];
rz(1.269066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3620944) q[2];
sx q[2];
rz(-2.4694314) q[2];
sx q[2];
rz(-0.54074311) q[2];
rz(2.8703459) q[3];
sx q[3];
rz(-1.2416779) q[3];
sx q[3];
rz(-1.1312243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75854492) q[0];
sx q[0];
rz(-0.31065148) q[0];
sx q[0];
rz(0.3824105) q[0];
rz(0.27579871) q[1];
sx q[1];
rz(-0.47960061) q[1];
sx q[1];
rz(-0.84024215) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58038515) q[0];
sx q[0];
rz(-1.2371644) q[0];
sx q[0];
rz(-1.946627) q[0];
x q[1];
rz(0.98767583) q[2];
sx q[2];
rz(-0.35735574) q[2];
sx q[2];
rz(0.59229702) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1237331) q[1];
sx q[1];
rz(-2.2869733) q[1];
sx q[1];
rz(-0.86114984) q[1];
rz(-1.6333079) q[3];
sx q[3];
rz(-2.0558236) q[3];
sx q[3];
rz(1.6597019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8648839) q[2];
sx q[2];
rz(-1.0397006) q[2];
sx q[2];
rz(-0.060982171) q[2];
rz(-1.7766772) q[3];
sx q[3];
rz(-2.958332) q[3];
sx q[3];
rz(2.0257559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0404496) q[0];
sx q[0];
rz(-1.018486) q[0];
sx q[0];
rz(-2.2818991) q[0];
rz(-1.0243833) q[1];
sx q[1];
rz(-2.5442217) q[1];
sx q[1];
rz(2.375367) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0939056) q[0];
sx q[0];
rz(-1.6062459) q[0];
sx q[0];
rz(0.043092502) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0417074) q[2];
sx q[2];
rz(-1.4204362) q[2];
sx q[2];
rz(-0.1783812) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.41118907) q[1];
sx q[1];
rz(-2.0918188) q[1];
sx q[1];
rz(-1.26544) q[1];
x q[2];
rz(-1.6556103) q[3];
sx q[3];
rz(-0.70422164) q[3];
sx q[3];
rz(-1.6617642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.031294558) q[2];
sx q[2];
rz(-0.86250192) q[2];
sx q[2];
rz(0.026570126) q[2];
rz(2.8734112) q[3];
sx q[3];
rz(-0.53332204) q[3];
sx q[3];
rz(-2.4466483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37102315) q[0];
sx q[0];
rz(-2.4257648) q[0];
sx q[0];
rz(0.41719607) q[0];
rz(-2.5542651) q[1];
sx q[1];
rz(-2.6081577) q[1];
sx q[1];
rz(0.74648285) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74325753) q[0];
sx q[0];
rz(-1.8881838) q[0];
sx q[0];
rz(2.7733735) q[0];
x q[1];
rz(2.3488147) q[2];
sx q[2];
rz(-0.93610349) q[2];
sx q[2];
rz(1.9605876) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.44297781) q[1];
sx q[1];
rz(-1.9483951) q[1];
sx q[1];
rz(2.5966132) q[1];
x q[2];
rz(1.0899441) q[3];
sx q[3];
rz(-1.1331416) q[3];
sx q[3];
rz(2.5646427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8314961) q[2];
sx q[2];
rz(-0.3336755) q[2];
sx q[2];
rz(-0.92158544) q[2];
rz(1.0725675) q[3];
sx q[3];
rz(-1.8662063) q[3];
sx q[3];
rz(0.5630365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4614586) q[0];
sx q[0];
rz(-2.7549094) q[0];
sx q[0];
rz(2.4899546) q[0];
rz(-0.98546511) q[1];
sx q[1];
rz(-0.54254222) q[1];
sx q[1];
rz(-2.998897) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089316018) q[0];
sx q[0];
rz(-2.8601649) q[0];
sx q[0];
rz(-1.3427585) q[0];
rz(-1.5161523) q[2];
sx q[2];
rz(-1.3846591) q[2];
sx q[2];
rz(-0.7344377) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.922077) q[1];
sx q[1];
rz(-2.0425046) q[1];
sx q[1];
rz(1.4336777) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0246932) q[3];
sx q[3];
rz(-1.6936692) q[3];
sx q[3];
rz(-1.6116796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5695213) q[2];
sx q[2];
rz(-0.60939747) q[2];
sx q[2];
rz(2.209668) q[2];
rz(-2.7052687) q[3];
sx q[3];
rz(-0.47567979) q[3];
sx q[3];
rz(1.9816192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6465004) q[0];
sx q[0];
rz(-2.2252872) q[0];
sx q[0];
rz(1.6012993) q[0];
rz(1.0013927) q[1];
sx q[1];
rz(-1.5512369) q[1];
sx q[1];
rz(-2.1839949) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0315899) q[0];
sx q[0];
rz(-2.6837807) q[0];
sx q[0];
rz(0.97719595) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6675435) q[2];
sx q[2];
rz(-1.821234) q[2];
sx q[2];
rz(-1.7360404) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5909429) q[1];
sx q[1];
rz(-2.6663482) q[1];
sx q[1];
rz(-1.5556504) q[1];
rz(-pi) q[2];
x q[2];
rz(0.05332975) q[3];
sx q[3];
rz(-0.55906534) q[3];
sx q[3];
rz(-0.4658162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6340948) q[2];
sx q[2];
rz(-1.8562506) q[2];
sx q[2];
rz(1.1320587) q[2];
rz(-3.0105528) q[3];
sx q[3];
rz(-1.9877501) q[3];
sx q[3];
rz(-1.8989782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89192724) q[0];
sx q[0];
rz(-0.51814336) q[0];
sx q[0];
rz(0.074935496) q[0];
rz(-0.19649188) q[1];
sx q[1];
rz(-2.7939929) q[1];
sx q[1];
rz(-0.67952716) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0655812) q[0];
sx q[0];
rz(-2.5437194) q[0];
sx q[0];
rz(-2.7779816) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6397298) q[2];
sx q[2];
rz(-0.40483311) q[2];
sx q[2];
rz(-0.73672494) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4008182) q[1];
sx q[1];
rz(-1.8535103) q[1];
sx q[1];
rz(-1.1307471) q[1];
rz(1.5854168) q[3];
sx q[3];
rz(-2.0818266) q[3];
sx q[3];
rz(2.8255812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3139265) q[2];
sx q[2];
rz(-2.1570666) q[2];
sx q[2];
rz(2.0619681) q[2];
rz(-0.45352724) q[3];
sx q[3];
rz(-0.87117666) q[3];
sx q[3];
rz(2.1760904) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6249348) q[0];
sx q[0];
rz(-0.17913945) q[0];
sx q[0];
rz(-3.0058885) q[0];
rz(-2.9756359) q[1];
sx q[1];
rz(-2.143492) q[1];
sx q[1];
rz(-0.96806324) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0128202) q[0];
sx q[0];
rz(-2.6261683) q[0];
sx q[0];
rz(1.5933871) q[0];
rz(-pi) q[1];
rz(2.9918475) q[2];
sx q[2];
rz(-2.1080842) q[2];
sx q[2];
rz(-1.6884691) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1727454) q[1];
sx q[1];
rz(-1.0732393) q[1];
sx q[1];
rz(-1.8109591) q[1];
x q[2];
rz(1.8289205) q[3];
sx q[3];
rz(-1.2069697) q[3];
sx q[3];
rz(0.3240594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.82219899) q[2];
sx q[2];
rz(-0.92331702) q[2];
sx q[2];
rz(-2.609002) q[2];
rz(-0.79832625) q[3];
sx q[3];
rz(-0.39755487) q[3];
sx q[3];
rz(-3.047191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.3548729) q[0];
sx q[0];
rz(-2.4479471) q[0];
sx q[0];
rz(-0.68674809) q[0];
rz(1.9101539) q[1];
sx q[1];
rz(-1.310692) q[1];
sx q[1];
rz(3.0795857) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3162295) q[0];
sx q[0];
rz(-2.4333409) q[0];
sx q[0];
rz(3.0407719) q[0];
rz(-pi) q[1];
rz(-0.19066452) q[2];
sx q[2];
rz(-1.4989509) q[2];
sx q[2];
rz(-2.004625) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1126642) q[1];
sx q[1];
rz(-0.24953574) q[1];
sx q[1];
rz(2.1965532) q[1];
rz(-pi) q[2];
rz(0.047895821) q[3];
sx q[3];
rz(-0.29828276) q[3];
sx q[3];
rz(-3.0259231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9318781) q[2];
sx q[2];
rz(-2.176602) q[2];
sx q[2];
rz(0.17267257) q[2];
rz(0.73575819) q[3];
sx q[3];
rz(-2.5685205) q[3];
sx q[3];
rz(2.6652523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19685766) q[0];
sx q[0];
rz(-1.4516964) q[0];
sx q[0];
rz(-1.4781937) q[0];
rz(-2.2810777) q[1];
sx q[1];
rz(-1.9445226) q[1];
sx q[1];
rz(-1.3938211) q[1];
rz(-0.71828193) q[2];
sx q[2];
rz(-0.16213308) q[2];
sx q[2];
rz(-2.8181974) q[2];
rz(-0.53027897) q[3];
sx q[3];
rz(-0.9835296) q[3];
sx q[3];
rz(0.26795878) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
