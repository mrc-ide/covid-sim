#include "Model.h"
#include <cstdlib>

bool Person::is_alive() const {
	return !this->is_dead();
}

bool Person::is_dead() const {
	return std::abs(this->inf) == InfStat_Dead;
}

void Person::make_susceptible() {
	this->inf = InfStat_Susceptible;
}

void Person::make_immune() {
	this->inf = InfStat_ImmuneAtStart;
}

void Person::make_infected() {
	this->inf = InfStat_Latent;
}

void Person::make_infectious_almost_symptomatic() {
	this->inf = InfStat_InfectiousAlmostSymptomatic;
}

void Person::make_infectious_asymptomatic() {
	this->inf = InfStat_InfectiousAsymptomaticNotCase;
}

void Person::make_case() {
	this->inf = InfStat_Case;
}

void Person::make_dead() {
	this->inf = (InfStat)(InfStat_Dead * this->inf / abs(this->inf));
}

void Person::make_recovered() {
	this->inf = (InfStat)(InfStat_Recovered * this->inf / std::abs(this->inf));
}

bool Person::is_susceptible() const {
	return this->inf == InfStat_Susceptible;
}

bool Person::is_latent() const {
	return this->inf == InfStat_Latent;
}

bool Person::is_symptomatic() const {
	return this->inf < InfStat_Susceptible;
}

bool Person::is_asymptomatic() const {
	return this->inf >= InfStat_Susceptible;
}

bool Person::is_case() const {
	return this->inf == InfStat_Case;
}

bool Person::is_not_case() const {
	return std::abs(this->inf) < InfStat_InfectiousAsymptomaticNotCase;
}

bool Person::is_infectious() const {
	return this->inf == InfStat_InfectiousAsymptomaticNotCase
		|| this->inf == InfStat_Case;
}

bool Person::is_infected() const {
	return this->inf != InfStat_Susceptible
		&& this->inf != InfStat_ImmuneAtStart;
}

bool Person::has_recovered() const {
	return std::abs(this->inf) == InfStat_Recovered;
}

bool Person::has_not_recovered() const {
	return std::abs(this->inf) < InfStat_Recovered;
}

bool Person::is_positive() const {
	return std::abs(this->inf) == InfStat_Case
		|| this->inf == InfStat_InfectiousAlmostSymptomatic;
}