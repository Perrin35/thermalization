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
rz(-1.0531955) q[0];
rz(1.6246417) q[1];
sx q[1];
rz(2.8487974) q[1];
sx q[1];
rz(8.4881633) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7335986) q[0];
sx q[0];
rz(-2.6331365) q[0];
sx q[0];
rz(1.5336799) q[0];
rz(-pi) q[1];
rz(2.1457315) q[2];
sx q[2];
rz(-0.81840063) q[2];
sx q[2];
rz(2.7374597) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9642311) q[1];
sx q[1];
rz(-1.291591) q[1];
sx q[1];
rz(-2.5482168) q[1];
rz(0.39150146) q[3];
sx q[3];
rz(-1.748852) q[3];
sx q[3];
rz(-0.72339971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8006353) q[2];
sx q[2];
rz(-0.95658797) q[2];
sx q[2];
rz(-1.0613649) q[2];
rz(-2.495885) q[3];
sx q[3];
rz(-2.784745) q[3];
sx q[3];
rz(2.1319353) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89479947) q[0];
sx q[0];
rz(-2.6428887) q[0];
sx q[0];
rz(2.6820768) q[0];
rz(-1.3872321) q[1];
sx q[1];
rz(-0.51001716) q[1];
sx q[1];
rz(-1.0646819) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1275528) q[0];
sx q[0];
rz(-1.8888374) q[0];
sx q[0];
rz(-1.7377338) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5881203) q[2];
sx q[2];
rz(-1.8380062) q[2];
sx q[2];
rz(1.653275) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.83512703) q[1];
sx q[1];
rz(-1.4773932) q[1];
sx q[1];
rz(-2.8207835) q[1];
rz(-pi) q[2];
rz(0.27402409) q[3];
sx q[3];
rz(-1.245387) q[3];
sx q[3];
rz(-2.7502378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.77949828) q[2];
sx q[2];
rz(-2.4694314) q[2];
sx q[2];
rz(2.6008495) q[2];
rz(0.27124673) q[3];
sx q[3];
rz(-1.2416779) q[3];
sx q[3];
rz(-2.0103683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.3830477) q[0];
sx q[0];
rz(-2.8309412) q[0];
sx q[0];
rz(-2.7591822) q[0];
rz(0.27579871) q[1];
sx q[1];
rz(-2.661992) q[1];
sx q[1];
rz(0.84024215) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58038515) q[0];
sx q[0];
rz(-1.9044283) q[0];
sx q[0];
rz(-1.946627) q[0];
rz(-pi) q[1];
rz(-1.2686549) q[2];
sx q[2];
rz(-1.7646175) q[2];
sx q[2];
rz(-0.42497488) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0178595) q[1];
sx q[1];
rz(-2.2869733) q[1];
sx q[1];
rz(2.2804428) q[1];
rz(1.5082848) q[3];
sx q[3];
rz(-2.0558236) q[3];
sx q[3];
rz(-1.4818907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.27670878) q[2];
sx q[2];
rz(-1.0397006) q[2];
sx q[2];
rz(3.0806105) q[2];
rz(-1.7766772) q[3];
sx q[3];
rz(-2.958332) q[3];
sx q[3];
rz(2.0257559) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0404496) q[0];
sx q[0];
rz(-2.1231066) q[0];
sx q[0];
rz(2.2818991) q[0];
rz(-2.1172093) q[1];
sx q[1];
rz(-2.5442217) q[1];
sx q[1];
rz(0.76622564) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1648772) q[0];
sx q[0];
rz(-3.0857997) q[0];
sx q[0];
rz(-2.4528422) q[0];
x q[1];
rz(-1.8623975) q[2];
sx q[2];
rz(-2.5935141) q[2];
sx q[2];
rz(-1.4983791) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.97570455) q[1];
sx q[1];
rz(-0.59671003) q[1];
sx q[1];
rz(2.6590682) q[1];
rz(-pi) q[2];
x q[2];
rz(0.071841876) q[3];
sx q[3];
rz(-2.2719682) q[3];
sx q[3];
rz(1.5506684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.031294558) q[2];
sx q[2];
rz(-2.2790907) q[2];
sx q[2];
rz(0.026570126) q[2];
rz(-2.8734112) q[3];
sx q[3];
rz(-0.53332204) q[3];
sx q[3];
rz(-0.69494438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-2.7705695) q[0];
sx q[0];
rz(-0.71582782) q[0];
sx q[0];
rz(-0.41719607) q[0];
rz(-0.5873276) q[1];
sx q[1];
rz(-0.53343499) q[1];
sx q[1];
rz(0.74648285) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3983351) q[0];
sx q[0];
rz(-1.8881838) q[0];
sx q[0];
rz(-2.7733735) q[0];
x q[1];
rz(-2.339614) q[2];
sx q[2];
rz(-2.1714513) q[2];
sx q[2];
rz(3.0026312) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.44297781) q[1];
sx q[1];
rz(-1.1931975) q[1];
sx q[1];
rz(-2.5966132) q[1];
x q[2];
rz(-2.0516486) q[3];
sx q[3];
rz(-1.1331416) q[3];
sx q[3];
rz(2.5646427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8314961) q[2];
sx q[2];
rz(-2.8079171) q[2];
sx q[2];
rz(-2.2200072) q[2];
rz(1.0725675) q[3];
sx q[3];
rz(-1.2753863) q[3];
sx q[3];
rz(-0.5630365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4614586) q[0];
sx q[0];
rz(-2.7549094) q[0];
sx q[0];
rz(2.4899546) q[0];
rz(0.98546511) q[1];
sx q[1];
rz(-0.54254222) q[1];
sx q[1];
rz(2.998897) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0522766) q[0];
sx q[0];
rz(-2.8601649) q[0];
sx q[0];
rz(1.3427585) q[0];
rz(-1.5161523) q[2];
sx q[2];
rz(-1.3846591) q[2];
sx q[2];
rz(-0.7344377) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.075293024) q[1];
sx q[1];
rz(-2.651804) q[1];
sx q[1];
rz(2.879786) q[1];
rz(-3.0050395) q[3];
sx q[3];
rz(-2.0210183) q[3];
sx q[3];
rz(3.0409851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5720713) q[2];
sx q[2];
rz(-0.60939747) q[2];
sx q[2];
rz(-2.209668) q[2];
rz(-2.7052687) q[3];
sx q[3];
rz(-0.47567979) q[3];
sx q[3];
rz(1.9816192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-0.49509224) q[0];
sx q[0];
rz(-0.91630542) q[0];
sx q[0];
rz(1.5402933) q[0];
rz(-2.1401999) q[1];
sx q[1];
rz(-1.5512369) q[1];
sx q[1];
rz(0.95759773) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1100028) q[0];
sx q[0];
rz(-0.45781198) q[0];
sx q[0];
rz(0.97719595) q[0];
rz(0.36105902) q[2];
sx q[2];
rz(-0.26810899) q[2];
sx q[2];
rz(-1.0323056) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5739096) q[1];
sx q[1];
rz(-2.0459818) q[1];
sx q[1];
rz(3.1337993) q[1];
rz(-pi) q[2];
x q[2];
rz(0.05332975) q[3];
sx q[3];
rz(-0.55906534) q[3];
sx q[3];
rz(2.6757765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5074978) q[2];
sx q[2];
rz(-1.2853421) q[2];
sx q[2];
rz(-2.0095339) q[2];
rz(-3.0105528) q[3];
sx q[3];
rz(-1.9877501) q[3];
sx q[3];
rz(-1.8989782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.2496654) q[0];
sx q[0];
rz(-2.6234493) q[0];
sx q[0];
rz(-3.0666572) q[0];
rz(0.19649188) q[1];
sx q[1];
rz(-0.34759977) q[1];
sx q[1];
rz(2.4620655) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20048902) q[0];
sx q[0];
rz(-1.3692432) q[0];
sx q[0];
rz(2.5747712) q[0];
x q[1];
rz(1.5018628) q[2];
sx q[2];
rz(-2.7367595) q[2];
sx q[2];
rz(0.73672494) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4008182) q[1];
sx q[1];
rz(-1.8535103) q[1];
sx q[1];
rz(1.1307471) q[1];
x q[2];
rz(1.5561759) q[3];
sx q[3];
rz(-1.059766) q[3];
sx q[3];
rz(-0.31601147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8276662) q[2];
sx q[2];
rz(-2.1570666) q[2];
sx q[2];
rz(2.0619681) q[2];
rz(-0.45352724) q[3];
sx q[3];
rz(-2.270416) q[3];
sx q[3];
rz(0.96550226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.6249348) q[0];
sx q[0];
rz(-0.17913945) q[0];
sx q[0];
rz(3.0058885) q[0];
rz(0.16595674) q[1];
sx q[1];
rz(-2.143492) q[1];
sx q[1];
rz(-0.96806324) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98685778) q[0];
sx q[0];
rz(-1.0555165) q[0];
sx q[0];
rz(-0.012795916) q[0];
rz(1.8161681) q[2];
sx q[2];
rz(-2.5858063) q[2];
sx q[2];
rz(-1.9751369) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.28169808) q[1];
sx q[1];
rz(-1.7813695) q[1];
sx q[1];
rz(2.631726) q[1];
rz(-0.59057136) q[3];
sx q[3];
rz(-0.44274032) q[3];
sx q[3];
rz(-0.31424943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3193937) q[2];
sx q[2];
rz(-2.2182756) q[2];
sx q[2];
rz(-2.609002) q[2];
rz(-0.79832625) q[3];
sx q[3];
rz(-0.39755487) q[3];
sx q[3];
rz(0.09440162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3548729) q[0];
sx q[0];
rz(-2.4479471) q[0];
sx q[0];
rz(0.68674809) q[0];
rz(-1.9101539) q[1];
sx q[1];
rz(-1.8309007) q[1];
sx q[1];
rz(-0.062006921) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3194763) q[0];
sx q[0];
rz(-1.6363167) q[0];
sx q[0];
rz(-2.4358552) q[0];
rz(-pi) q[1];
rz(1.4976296) q[2];
sx q[2];
rz(-1.760963) q[2];
sx q[2];
rz(-0.41997465) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0289284) q[1];
sx q[1];
rz(-0.24953574) q[1];
sx q[1];
rz(2.1965532) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5560772) q[3];
sx q[3];
rz(-1.2728661) q[3];
sx q[3];
rz(3.0760279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9318781) q[2];
sx q[2];
rz(-2.176602) q[2];
sx q[2];
rz(0.17267257) q[2];
rz(-0.73575819) q[3];
sx q[3];
rz(-0.57307214) q[3];
sx q[3];
rz(-0.47634038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19685766) q[0];
sx q[0];
rz(-1.6898962) q[0];
sx q[0];
rz(1.663399) q[0];
rz(-0.860515) q[1];
sx q[1];
rz(-1.1970701) q[1];
sx q[1];
rz(1.7477716) q[1];
rz(2.4233107) q[2];
sx q[2];
rz(-0.16213308) q[2];
sx q[2];
rz(-2.8181974) q[2];
rz(0.53027897) q[3];
sx q[3];
rz(-2.1580631) q[3];
sx q[3];
rz(-2.8736339) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
