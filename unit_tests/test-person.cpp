#include <gtest/gtest.h>

#include "Models/Person.h"

TEST(CovidSimPersonTests, Person_Susceptible)
{
  Person* P = new Person();
  
  P->set_susceptible();
  ASSERT_FALSE(P->do_not_vaccinate());
  ASSERT_FALSE(P->is_dead());
  ASSERT_FALSE(P->is_dead_was_asymp());
  ASSERT_FALSE(P->is_dead_was_symp());
  ASSERT_FALSE(P->is_case());
  ASSERT_FALSE(P->is_immune_at_start());
  ASSERT_FALSE(P->is_infectious_almost_symptomatic());
  ASSERT_FALSE(P->is_infectious_asymptomatic_not_case());
  ASSERT_FALSE(P->is_latent());
  ASSERT_FALSE(P->is_never_symptomatic());
  ASSERT_FALSE(P->is_recovered());
  ASSERT_FALSE(P->is_recovered_symp());
  ASSERT_TRUE(P->is_alive());
  ASSERT_TRUE(P->is_not_yet_symptomatic());
  ASSERT_TRUE(P->is_susceptible_or_infected());
  delete P;
}

TEST(CovidSimPersonTests, Person_Latent)
{
  Person* P = new Person();
  P->set_latent();
  ASSERT_FALSE(P->do_not_vaccinate());
  ASSERT_FALSE(P->is_dead());
  ASSERT_FALSE(P->is_dead_was_asymp());
  ASSERT_FALSE(P->is_dead_was_symp());
  ASSERT_FALSE(P->is_case());
  ASSERT_FALSE(P->is_immune_at_start());
  ASSERT_FALSE(P->is_infectious_almost_symptomatic());
  ASSERT_FALSE(P->is_infectious_asymptomatic_not_case());
  ASSERT_TRUE(P->is_latent());
  ASSERT_TRUE(P->is_never_symptomatic());
  ASSERT_FALSE(P->is_recovered());
  ASSERT_FALSE(P->is_recovered_symp());
  ASSERT_TRUE(P->is_alive());
  ASSERT_TRUE(P->is_not_yet_symptomatic());
  ASSERT_TRUE(P->is_susceptible_or_infected());

  delete P;
}

TEST(CovidSimPersonTests, Person_Inf_Alm_Sympt)
{
  Person* P = new Person();
  P->set_infectious_almost_symptomatic();
  ASSERT_FALSE(P->do_not_vaccinate());
  ASSERT_FALSE(P->is_dead());
  ASSERT_FALSE(P->is_dead_was_asymp());
  ASSERT_FALSE(P->is_dead_was_symp());
  ASSERT_FALSE(P->is_case());
  ASSERT_FALSE(P->is_immune_at_start());
  ASSERT_TRUE(P->is_infectious_almost_symptomatic());
  ASSERT_FALSE(P->is_infectious_asymptomatic_not_case());
  ASSERT_FALSE(P->is_latent());
  ASSERT_FALSE(P->is_never_symptomatic());
  ASSERT_FALSE(P->is_recovered());
  ASSERT_FALSE(P->is_recovered_symp());
  ASSERT_TRUE(P->is_alive());
  ASSERT_TRUE(P->is_not_yet_symptomatic());
  ASSERT_TRUE(P->is_susceptible_or_infected());
  delete P;
}

TEST(CovidSimPersonTests, Person_Case)
{
  Person* P = new Person();
  P->set_case();
  ASSERT_TRUE(P->do_not_vaccinate());
  ASSERT_FALSE(P->is_dead());
  ASSERT_FALSE(P->is_dead_was_asymp());
  ASSERT_FALSE(P->is_dead_was_symp());
  ASSERT_TRUE(P->is_case());
  ASSERT_FALSE(P->is_immune_at_start());
  ASSERT_FALSE(P->is_infectious_almost_symptomatic());
  ASSERT_FALSE(P->is_infectious_asymptomatic_not_case());
  ASSERT_FALSE(P->is_latent());
  ASSERT_FALSE(P->is_never_symptomatic());
  ASSERT_FALSE(P->is_recovered());
  ASSERT_FALSE(P->is_recovered_symp());
  ASSERT_TRUE(P->is_alive());
  ASSERT_FALSE(P->is_not_yet_symptomatic());
  ASSERT_TRUE(P->is_susceptible_or_infected());
  delete P;
}

TEST(CovidSimPersonTests, Person_Inf_Asymp_Not_Case)
{
  Person* P = new Person();
  P->set_infectious_asymptomatic_not_case();
  ASSERT_FALSE(P->do_not_vaccinate());
  ASSERT_FALSE(P->is_dead());
  ASSERT_FALSE(P->is_dead_was_asymp());
  ASSERT_FALSE(P->is_dead_was_symp());
  ASSERT_FALSE(P->is_case());
  ASSERT_FALSE(P->is_immune_at_start());
  ASSERT_FALSE(P->is_infectious_almost_symptomatic());
  ASSERT_TRUE(P->is_infectious_asymptomatic_not_case());
  ASSERT_FALSE(P->is_latent());
  ASSERT_TRUE(P->is_never_symptomatic());
  ASSERT_FALSE(P->is_recovered());
  ASSERT_FALSE(P->is_recovered_symp());
  ASSERT_TRUE(P->is_alive());
  ASSERT_FALSE(P->is_not_yet_symptomatic());
  ASSERT_TRUE(P->is_susceptible_or_infected());
  delete P;
}

TEST(CovidSimPersonTests, Person_Death_Symp)
{
  Person* P = new Person();
  P->set_case();
  P->set_dead();
  ASSERT_TRUE(P->do_not_vaccinate());
  ASSERT_TRUE(P->is_dead());
  ASSERT_FALSE(P->is_dead_was_asymp());
  ASSERT_TRUE(P->is_dead_was_symp());
  ASSERT_FALSE(P->is_case());
  ASSERT_FALSE(P->is_immune_at_start());
  ASSERT_FALSE(P->is_infectious_almost_symptomatic());
  ASSERT_FALSE(P->is_infectious_asymptomatic_not_case());
  ASSERT_FALSE(P->is_latent());
  ASSERT_FALSE(P->is_never_symptomatic());
  ASSERT_FALSE(P->is_recovered());
  ASSERT_FALSE(P->is_recovered_symp());
  ASSERT_FALSE(P->is_alive());
  ASSERT_FALSE(P->is_not_yet_symptomatic());
  ASSERT_FALSE(P->is_susceptible_or_infected());
  delete P;
}

TEST(CovidSimPersonTests, Person_Death_ASymp)
{
  Person* P = new Person();
  P->set_susceptible();
  P->set_dead();
  ASSERT_TRUE(P->do_not_vaccinate());
  ASSERT_TRUE(P->is_dead());
  ASSERT_TRUE(P->is_dead_was_asymp());
  ASSERT_FALSE(P->is_dead_was_symp());
  ASSERT_FALSE(P->is_case());
  ASSERT_FALSE(P->is_immune_at_start());
  ASSERT_FALSE(P->is_infectious_almost_symptomatic());
  ASSERT_FALSE(P->is_infectious_asymptomatic_not_case());
  ASSERT_FALSE(P->is_latent());
  ASSERT_TRUE(P->is_never_symptomatic());
  ASSERT_FALSE(P->is_recovered());
  ASSERT_FALSE(P->is_recovered_symp());
  ASSERT_FALSE(P->is_alive());
  ASSERT_FALSE(P->is_not_yet_symptomatic());
  ASSERT_FALSE(P->is_susceptible_or_infected());
  delete P;
}

TEST(CovidSimPersonTests, Person_Recover_Symp)
{
  Person* P = new Person();
  P->set_case();
  P->set_recovered();
  ASSERT_TRUE(P->do_not_vaccinate());
  ASSERT_FALSE(P->is_dead());
  ASSERT_FALSE(P->is_dead_was_asymp());
  ASSERT_FALSE(P->is_dead_was_symp());
  ASSERT_FALSE(P->is_case());
  ASSERT_FALSE(P->is_immune_at_start());
  ASSERT_FALSE(P->is_infectious_almost_symptomatic());
  ASSERT_FALSE(P->is_infectious_asymptomatic_not_case());
  ASSERT_FALSE(P->is_latent());
  ASSERT_FALSE(P->is_never_symptomatic());
  ASSERT_TRUE(P->is_recovered());
  ASSERT_TRUE(P->is_recovered_symp());
  ASSERT_TRUE(P->is_alive());
  ASSERT_FALSE(P->is_not_yet_symptomatic());
  ASSERT_FALSE(P->is_susceptible_or_infected());
  delete P;
}

TEST(CovidSimPersonTests, Person_Recover_ASymp)
{
  Person* P = new Person();
  P->set_susceptible();
  P->set_recovered();
  ASSERT_FALSE(P->do_not_vaccinate());
  ASSERT_FALSE(P->is_dead());
  ASSERT_FALSE(P->is_dead_was_asymp());
  ASSERT_FALSE(P->is_dead_was_symp());
  ASSERT_FALSE(P->is_case());
  ASSERT_FALSE(P->is_immune_at_start());
  ASSERT_FALSE(P->is_infectious_almost_symptomatic());
  ASSERT_FALSE(P->is_infectious_asymptomatic_not_case());
  ASSERT_FALSE(P->is_latent());
  ASSERT_TRUE(P->is_never_symptomatic());
  ASSERT_TRUE(P->is_recovered());
  ASSERT_FALSE(P->is_recovered_symp());
  ASSERT_TRUE(P->is_alive());
  ASSERT_FALSE(P->is_not_yet_symptomatic());
  ASSERT_FALSE(P->is_susceptible_or_infected());
  delete P;
}
